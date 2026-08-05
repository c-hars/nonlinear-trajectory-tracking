function [u,solve_info] = compute_u_SDDRE_v3(tk,xk,k,uk,r_,C,Qy,R,Qyf,tf,qp,opts)
    arguments
        tk
        xk
        k
        uk
        r_
        C
        Qy
        R
        Qyf
        tf
        qp
        opts.DecimationFactor = 1 % idea - for future work :) Not implemented.
        opts.SDC_A_function = @ get_A_matrix_SDRE_EulerAttitude
        opts.SDC_B_function = @ get_B_matrix_SDRE
        opts.PreviewHorizon = 2.0 % this is all really designed for finite-horizon formulations, with the horizon large (>= 2.0 for our case) -- but inf is a safe choice, just get delayed tracking
        opts.UseFullFiniteHorizonMPCAtTerminal = true % best set to false if you need precisely consistent/predictable solve times - adds a bit of overhead to do the full recursion for K too
        opts.AlwaysUseFullFiniteHorizonMPC = false % set to true → regression to the standard MPC cost function being optimised / directly comparable OCP with LinMPC, NLMPC -- no longer optimal preview control that uses the infinite/finite horizon split (core part of the algorithm!)

        opts.DARESolver (1,:) char {mustBeMember(opts.DARESolver, {'nk','riccati','cold','idare','dlqr','sda'})} = 'nk'
        opts.DARESolverOpts (1,1) struct = struct()
    end

    % --- 1. Compute the (ZOH-discretised) SDC matrices ---

    clock_start = tic;
    Ac = opts.SDC_A_function(xk, qp);
    uk = max(uk, -qp.nominal_omegas(:));
    uk = min(uk,  qp.max_du(:));
    Bc = opts.SDC_B_function(uk, qp);
    [A,B] = c2d_zoh_expm(Ac,Bc,qp.Ts);
    solve_info.TimeToComputeDiscreteSDCMatrices = toc(clock_start);


    % --- 2. Solve the DARE for P_ss, K_ss ---

    persistent P_ss K_ss
    persistent warnedPreviously % for loud fallback warnings
    if k == 1, warnedPreviously = false; end

    clock_start = tic;

    if k == 1, opts.DARESolver = 'cold'; end  % cold solve on init

    switch lower(opts.DARESolver)

    case {'cold','sda'}
        sdaOpts = struct();
        if isfield(opts.DARESolverOpts, 'Tolerance')
            sdaOpts.Tolerance = opts.DARESolverOpts.Tolerance;
        end
        args = namedargs2cell(sdaOpts);
        [P_ss, info] = dare_sda(A, B, C'*Qy*C, R, args{:});
        K_ss = (R + B'*P_ss*B) \ (B'*P_ss*A);

        solve_info.DARESolverTolAchieved = info.TolAchieved;
        solve_info.DARESolverNumIters    = info.SolverIterations;
        solve_info.DARESolverSuccess     = info.SolveSuccess;

    case 'idare'

        [P_ss,K_ss,~,info]  = idare(A, B, C'*Qy*C, R);
        parse_idare_info(info);

        solve_info.DARESolverTolAchieved = compute_dare_residual(A,B,C'*Qy*C,R,P_ss,K_ss);
        solve_info.DARESolverNumIters    = 0;
        solve_info.DARESolverSuccess     = 1;

    case 'dlqr'

        [K_ss,P_ss] = dlqr(A, B, C'*Qy*C, R);
        solve_info.DARESolverTolAchieved = compute_dare_residual(A,B,C'*Qy*C,R,P_ss,K_ss);
        solve_info.DARESolverNumIters    = 0;
        solve_info.DARESolverSuccess     = 1;

    case {'nk','riccati'} % the two iterative methods

        solverOpts = opts.DARESolverOpts;
        solverOpts.Method = opts.DARESolver;
        args = namedargs2cell(solverOpts);
        [P_ss, info] = iterative_dare(A, B, C'*Qy*C, R, P_ss, args{:});
        K_ss  = (R + B'*P_ss*B) \ (B'*P_ss*A);

        % Fallback to SDA if initialised outside the stability basin (only needed for NK)
        solve_info.DARESolverSuccess = info.SolveSuccess;
        if strcmpi(opts.DARESolver,'nk') && ~info.SolveSuccess
            if isfield(info,'UnstableK0') && info.UnstableK0

                sdaOpts = struct();
                if isfield(opts.DARESolverOpts, 'Tolerance'), sdaOpts.Tolerance = opts.DARESolverOpts.Tolerance; end
                args = namedargs2cell(sdaOpts);
                [P_ss, info] = dare_sda(A, B, C'*Qy*C, R, args{:});
                K_ss = (R + B'*P_ss*B) \ (B'*P_ss*A);
                
                solve_info.DARESolverSuccess = 0.5*info.SolveSuccess; % 0.5 == sentinel value for partial success
                
                if warnedPreviously == false
                    warning("compute_u_SDDRE_v3: NK iteration initialised outside of stability basin at k=%d. Fell back to SDA cold solve: success flag was %.1f", k, solve_info.DARESolverSuccess)
                    warnedPreviously = true;
                end
                
            end
        end

        solve_info.DARESolverTolAchieved = info.TolAchieved;
        solve_info.DARESolverNumIters    = info.SolverIterations;

    otherwise

        error("Unknown DARESolver type.");

    end
    solve_info.TimeToSolveDARE = toc(clock_start);
    

    % --- 3. Compute the feedforward (reference preview) term ---
    
    clock_start = tic;
    A_cl  = A - B*K_ss;
    if isinf(opts.PreviewHorizon)
        % Constant reference approximation
        s_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*s_k1);
    else
        M = round(opts.PreviewHorizon/qp.Ts);  % receding horizon preview window
        approachingTerminal = (k+M >= length(r_));
        
        if ~(opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal))
            % Main branch: infinite horizon w/ finite preview
        
            % v_k1 = C'*Qy*r_(:, min(k+M,length(r_))); % normal costate seed
            v_k1 = (eye(12) - A_cl') \ (C'*Qy*r_(:, min(k+M,length(r_))));
            for j = M-1:-1:1
                v_k1 = A_cl'*v_k1 + C'*Qy*r_(:, min(k+j,length(r_)));
            end
            u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*v_k1);
        else
            % Approaching terminal state.
            % Compute full recursion for K_k, K_k^v, v_{k+1}.
            % This ensures K is ramped up appropriately according to Qyf.
            % (This is then exactly equivalent to LQT MPC using frozen A_SDC(xk), B_SDC(xk))
            P = C'*Qyf*C;
            M = min(M, length(r_) - k);
            v = C'*Qyf*r_(:, k+M);
            for j = M-1:-1:1
                K_j    = (R + B'*P*B) \ (B'*P*A);
                A_cl_j = A - B*K_j;
                P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                v      = A_cl_j'*v + C'*Qy*r_(:, k+j);
            end
            Kk  = (R + B'*P*B) \ (B'*P*A);
            Kvk = (R + B'*P*B) \ B';
            u   = -Kk*xk + Kvk*v;
        end
    end
    solve_info.TimeToComputeFeedforward = toc(clock_start);
end


function parse_idare_info(info)
    switch info.Report
        case 1, warning("idare(), info.Report == 1 (The solution accuracy is poor)")
        case 2, warning("idare(), info.Report == 2 (The solution is not finite)")
        case 3, error("idare(), info.Report == 3 (No solution found since the Symplectic spectrum, denoted by [L;1./L], has eigenvalues on the unit circle)")
        case 4, error("idare(), info.Report == 4 (Pencil is singular ([B;S;R] is rank deficient)")
    end
end