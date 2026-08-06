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

                warning("compute_u_SDDRE_v3: NK iteration initialised outside of stability basin at k=%d. Fell back to SDA cold solve: success flag was %.1f", k, solve_info.DARESolverSuccess)
                
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

    persistent CtQyr_
    % Weighted reference, CtQy*r_
    % This quantity is referenced throughout - for efficiency, prefer precomputing as a fixed offline table. Recomputing it each call is redundant computation.
    % But if the reference plan changes, or you're feeding a constantly changing reference, not a known-apriori plan - then you can't get away with precomputing. A circular buffer sized to the preview window (and only updated with the one-step new additions) would be a good idea there.
    % Padded by M+d so no index clamping is needed
    if k == 1 || isempty(CtQyr_)
        d = max(1, round(opts.DecimationFactor));
        M = round(opts.PreviewHorizon / qp.Ts);
        rpad = [r_, repmat(r_(:,end), 1, M + d)];
        CtQyr_ = C' * Qy * rpad;
    end

    if isinf(opts.PreviewHorizon)
        % Constant reference approximation
        s_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*s_k1);

    else
        M = round(opts.PreviewHorizon/qp.Ts);   % receding horizon preview window
        approachingTerminal = (k+M >= length(r_));
        d = max(1, round(opts.DecimationFactor));

        if ~(opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal))

            F = A_cl';

            if d == 1
                v_k1 = (eye(12) - F) \ CtQyr_(:, k+M);
                for j = M-1:-1:1
                    v_k1 = F*v_k1 + CtQyr_(:, k+j);
                end
            else

                Ginf = (eye(12) - F) \ eye(12); % sum_{i=0 to inf} F^i
                Fd = F^d;
                Gd   = (eye(12) - Fd) * Ginf;   % sum_{i=0 to d-1} F^i

                nsteps = M - 1;
                nearFieldFineSteps = 1*d; % N blocks, always at full res
                % nearFieldFineSteps = ceil(0.05/qp.Ts); % 0.05s, always at full res
                nfine  = min(nsteps, nearFieldFineSteps);
                nblk   = floor((nsteps - nfine)/d);
                nexact = nsteps - nblk*d;  % >= nfine by construction

                v_k1 = (eye(12) - F) \ CtQyr_(:, k+M);
                % far field
                for b = nblk:-1:1
                    jblk = nexact + (b-1)*d + 1;  % the first step in block b
                    v_k1 = Fd*v_k1 + Gd*mean(CtQyr_(:, k+jblk : k+jblk+d-1), 2);
                end
                % near field remnant
                for j = nexact:-1:1
                    v_k1 = F*v_k1 + CtQyr_(:, k+j);
                end
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

            if d == 1

                for j = M-1:-1:1
                    K_j    = (R + B'*P*B) \ (B'*P*A);
                    A_cl_j = A - B*K_j;
                    P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                    v      = A_cl_j'*v + CtQyr_(:, k+j);
                end

            else

                % [Ad,Bd] = c2d_zoh_expm(Ac,Bc,qp.Ts * d);
                %   NB: Ad != "A discrete", it's "A decimated" (but yes, discrete time also)
                % Below uses a cheaper re-discretisation (using the already-computed matrix exponential).
                % > Exact ZOH discretisation at d*Ts, built from the fine-step (A,B)  
                % > Taking d fine steps with a constant input is the same as one coarse step
                %   with Ad = A^d and Bd = (I + A + A^2 + ... + A^{d-1}) * B = Gd * B
                Apow = eye(12);
                Gd   = eye(12);
                for i = 1:d-1
                    Apow = A * Apow;
                    Gd   = Gd + Apow;
                end
                Ad = A * Apow;
                Bd = Gd * B;

                Qyd = Qy*d;
                Rd  = R*d;

                nsteps = M-1;
                nearFieldFineSteps = 1*d; % N blocks, always at full res
                % nearFieldFineSteps = ceil(0.05/qp.Ts); % 0.05s, always at full res
                nfine  = min(nsteps, nearFieldFineSteps); % steps adjacent to k forced to full Ts
                nblk   = floor((nsteps - nfine)/d);
                nexact = nsteps - nblk*d;  % >= nfine by construction

                % "blocks" at coarser, decimated Ts
                for b = nblk:-1:1
                    jblk   = nexact + (b-1)*d + 1;  % the first step in block b
                    K_b    = (Rd + Bd'*P*Bd) \ (Bd'*P*Ad);
                    A_cl_b = Ad - Bd*K_b;
                    P      = C'*Qyd*C + K_b'*Rd*K_b + A_cl_b'*P*A_cl_b;
                    v = A_cl_b'*v + d*mean(CtQyr_(:, k+jblk : k+jblk+d-1), 2);
                end
                % near field remnant
                % same as other branch, done at the full Ts
                for j = nexact:-1:1
                    K_j    = (R + B'*P*B) \ (B'*P*A);
                    A_cl_j = A - B*K_j;
                    P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                    v      = A_cl_j'*v + CtQyr_(:, k+j);
                end

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