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
    CtQy  = C'*Qy;

    if isinf(opts.PreviewHorizon)
        % Constant reference approximation
        s_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*s_k1);

    else
        M = round(opts.PreviewHorizon/qp.Ts);   % receding horizon preview window
        approachingTerminal = (k+M >= length(r_));

        % Weighted reference, CtQy*r_
        % Rebuilt each call for clarity. In implementations you'd probably want to precompute this as a fixed offline table - recomputing it each call is redundant computation
        % Padded by M+d so no index clamping is needed
        d = max(1, round(opts.DecimationFactor));
        rpad = [r_, repmat(r_(:,end), 1, M + d)];
        CtQyr_  = CtQy * rpad;

        if ~(opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal))

            F = A_cl';

            if d == 1
                v_k1 = (eye(12) - F) \ CtQyr_(:, k+M);
                for j = M-1:-1:1
                    v_k1 = F*v_k1 + CtQyr_(:, k+j);
                end
            else

                Ginf = (eye(12) - F) \ eye(12);
                Fd = F^d;
                Gd   = (eye(12) - Fd) * Ginf;  % sum_{i<d} F^i

                nsteps = M - 1;
                nblk   = floor(nsteps/d);
                nexact = nsteps - nblk*d;

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
                    v      = A_cl_j'*v + C'*Qy*r_(:, k+j);
                end

            else

                % The below is not great, too analogous porting. Wouldve been fine if it worked but it doesnt anyays.
                % Cleaner: just discretise at the decimated Ts lol, run the recursion there. That's the whole point of the decimation factor

                max_dx = [[1 1 1]*10*0.01, ...      % | 10 | 5   | [cm]
                          [1 1 1]*5*0.01, ...       % | 5  | 5   | [cm/s]
                          deg2rad([1 1 1]*5), ...   % | 5  | 2.5 | [deg]
                          deg2rad([1 1 1]*50)]';    % | 50 | 50  | [deg/s]
                jstop = 1; % (*)
                epsstop = 0.2; % (*) delta u threshold
                checkEvery = 25; % (*) it's wasteful to check every iteration
                checkPast = 300;  % (*) also wasteful to check too early
                cnt = 0;
                for j = M-1:-1:1
                    K_j    = (R + B'*P*B) \ (B'*P*A);
                    A_cl_j = A - B*K_j;
                    P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                    v      = A_cl_j'*v + C'*Qy*r_(:, k+j);
                    if cnt > checkPast && mod(cnt,checkEvery) == 0 % (*)
                        if norm((K_j - K_prev)*max_dx,'fro') <= epsstop % (*)
                            jstop = j; % Riccati has settled; rest is constant-F % (*)
                            % fprintf("Settled after %d iters\n", M-1-j) % (*)
                            disp((j-(M-1))*-1) % (*)
                            break % (*)
                        end % (*)
                    end % (*)
                    K_prev = K_j; % (*)
                    cnt = cnt + 1; % (*)
                end

                nsteps = jstop - 1;  % (*)
                if nsteps > 0  % (*)
                    F    = (A - B*K_ss)';  % (*)
                    % F = A_cl_j';
                    Ginf = (eye(12) - F) \ eye(12);  % (*)
                    Fd   = F^d;  % (*)
                    Gd   = (eye(12) - Fd) * Ginf;       % sum_{i<d} F^i  % (*)

                    nblk   = floor(nsteps/d);  % (*)
                    nexact = nsteps - nblk*d;  % (*)

                    for b = nblk:-1:1  % (*)
                        jblk = nexact + (b-1)*d + 1;  % (*)
                        v = Fd*v + Gd*mean(CtQyr_(:, k+jblk : k+jblk+d-1), 2);  % (*)
                    end  % (*)
                    for j = nexact:-1:1  % (*)
                        v = F*v + CtQyr_(:, k+j);  % (*)
                    end  % (*)
                end  % (*)
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