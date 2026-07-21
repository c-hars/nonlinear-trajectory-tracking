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
        opts.PreviewHorizon = 2.0 % this is all really designed for finite-horizon formulations, with the horizon large (>= 2.0 for our case) -- but inf is a safe default choice, just get delayed tracking
        opts.UseFullFiniteHorizonMPCAtTerminal = true % best set to false if you need precisely consistent/predictable solve times - adds a bit of overhead to do the recursion for K too
        opts.AlwaysUseFullFiniteHorizonMPC = false % set to false → the feedback part, u=Kx, does not see the terminal condition, instead relying on the infinite horizon gain, maintained via P_ss, resolved via warm start C-NK

        opts.DARESolver (1,:) char {mustBeMember(opts.DARESolver, {'nk','riccati','cold','idare','dlqr'})} = 'nk'
        opts.DARESolverOpts (1,1) struct = struct()
    end

    clock_start = tic;
    A = opts.SDC_A_function(xk, qp);
    uk = max(uk, -qp.nominal_omegas(:));
    uk = min(uk,  qp.max_du(:));
    B = opts.SDC_B_function(uk, qp);
    [A,B] = c2d_zoh_expm(A,B,qp.Ts);
    solve_info.TimeToComputeDiscreteSDCMatrices = toc(clock_start);


    persistent P_ss K_ss
    persistent warnedPreviously % for loud fallback warnings

    clock_start = tic;
    if k == 1 || strcmp(opts.DARESolver, 'cold') || strcmp(opts.DARESolver, 'idare')

        [P_ss,K_ss,~,info]  = idare(A, B, C'*Qy*C, R);
        switch info.Report
            case 1, warning("idare(), info.Report == 1 (The solution accuracy is poor)")
            case 2, warning("idare(), info.Report == 2 (The solution is not finite)")
            case 3, disp(A); error("idare(), info.Report == 3 (No solution found since the Symplectic spectrum, denoted by [L;1./L], has eigenvalues on the unit circle)")
            case 4, error("idare(), info.Report == 4 (Pencil is singular ([B;S;R] is rank deficient)")
        end

        solve_info.DARESolverTolAchieved = compute_dare_residual(A,B,C'*Qy*C,R,P_ss,K_ss);
        solve_info.DARESolverNumIters    = 0;
        solve_info.DARESolverSuccess     = 1;

    elseif strcmpi(opts.DARESolver, 'dlqr')

        [K_ss,P_ss] = dlqr(A, B, C'*Qy*C, R);
        solve_info.DARESolverTolAchieved = compute_dare_residual(A,B,C'*Qy*C,R,P_ss,K_ss);
        solve_info.DARESolverNumIters    = 0;
        solve_info.DARESolverSuccess     = 1;

    else

        solverOpts = opts.DARESolverOpts;
        solverOpts.Method = opts.DARESolver;
        args = namedargs2cell(solverOpts);
        [P_ss, info] = iterative_dare(A, B, C'*Qy*C, R, P_ss, args{:});

        solve_info.DARESolverTolAchieved = info.TolAchieved;  % no need to compute_dare_residual() again
        solve_info.DARESolverNumIters    = info.SolverIterations;
        solve_info.DARESolverSuccess     = info.SolveSuccess;
        
        K_ss  = (R + B'*P_ss*B) \ (B'*P_ss*A);

        % This is a possibility with NK...
        % NK iteration returns a destabilising gain if K0 is destabilising.
        % For future work, an eigenvalue check before NK iteration - for now, fallback with loud warning.
        if ~all(abs(eig(A-B*K_ss)) < 1)
            if isempty(warnedPreviously) || warnedPreviously == false
                warning("k=%d: iterative dare returned destabilising gain, reverting to Riccati", k);
                warnedPreviously = true;
            end

            [P_ss, info] = iterative_dare(A, B, C'*Qy*C, R, P_ss, Method='Riccati', MinIters=25, MaxIters=25*100, Tolerance=1e-4);
    
            solve_info.DARESolverTolAchieved = info.TolAchieved;  % no need to compute_dare_residual() again
            solve_info.DARESolverNumIters    = info.SolverIterations;
            solve_info.DARESolverSuccess     = 0;

            K_ss  = (R + B'*P_ss*B) \ (B'*P_ss*A);
        end

    end
    solve_info.TimeToSolveDARE = toc(clock_start);
    

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