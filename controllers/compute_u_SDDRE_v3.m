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
        opts.DecimationFactors = 1 % scalar or ascending vector of decimation factors. Tiers 1..end-1 are near-horizon (each gets StepsPerTier coarse steps); last element is the far-horizon tier absorbing the rest.
        opts.StepsPerTier      = 1 % number of coarse steps per near-horizon tier
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
    % A = eye(12) + Ac*qp.Ts;
    % B = Bc*qp.Ts;
    solve_info.TimeToComputeDiscreteSDCMatrices = toc(clock_start);


    % --- 2. Solve the DARE for P_ss, K_ss ---

    clock_start = tic;
    Q = C'*Qy*C;
    [P_ss, K_ss, dare_info] = solve_dare(A, B, Q, R, k, opts.DARESolver, opts.DARESolverOpts);
    solve_info.DARESolverTolAchieved = dare_info.TolAchieved;
    solve_info.DARESolverNumIters    = dare_info.NumIters;
    solve_info.DARESolverSuccess     = dare_info.Success;
    solve_info.TimeToSolveDARE = toc(clock_start);
    

    % --- 3. Compute the feedforward (preview) term ---
    
    clock_start = tic;
    A_cl  = A - B*K_ss;

    if isinf(opts.PreviewHorizon)
        % Constant reference approximation
        v_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*v_k1);

    else
        n_preview = round(opts.PreviewHorizon/qp.Ts);   % receding horizon preview window
        approachingTerminal = (k + n_preview >= length(r_));


        % --- Preview grid terminology ---
        % 
        % Fine step: a prediction step of duration Ts.
        %
        % Coarse step: a prediction step covering more time than a fine step (duration > Ts). Lower model fidelity and accuracy, but faster compute (in aggregate). A coarse step has duration d*Ts, where `d` is the decimation factor; `d` consecutive fine steps are grouped and approximated as one coarse step.
        %
        % Block: each discrete-time step is over d*Ts seconds -- 'block' is shorthand for this chunk of time.
        % 
        % Coarse system matrices: (A_d,B_d,Qy_d,R_d) - the system matrices (A,B,Q,R) after applying a decimation factor. Represents the same underlying system dynamics (Ac,Bc) and OCP formulation (Qy_c,R_c), but discretised at a different (lower-resolution) sampling rate.
        %
        % Fine/coarse recursion: a backward recursion processing one fine/coarse step at a time.
        % 
        % Tier: a contiguous sequence of steps sharing the same decimation factor.
        % - Far-horizon tier: the final tier -- farthest from the current time, and uses the largest decimation factor. Absorbs whatever horizon remains after the near-horizon tiers.
        % - Near-horizon tier: any tier other than the far-horizon tier. Each holds StepsPerTier steps.
        %
        % Remainder fine steps: fine steps nearest to the current time that don't fit into the tier schedule; processed by the fine recursion.
        %
        % Preview grid: the full arrangement across the preview horizon. Forward in time:
        %    remainder fine steps -> near-horizon tiers -> far-horizon tier -> preview endpoint
        % The backward recursion processes this in reverse.



        % Preview grid setup. Shared by both branches.
        df_ = sort(opts.DecimationFactors(:)', 'ascend');
        n_coarse_per_tier = opts.StepsPerTier;
        n_tiers = length(df_);

        if ~(opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal))

            F = A_cl';
            n_fine = n_preview - 1;
            n_fine_near = n_coarse_per_tier * sum(df_(1:end-1));

            if n_fine_near >= n_fine
                % Fall back to fine recursion
                if k == 1
                    warning('compute_u_SDDRE_v3: near-horizon tiers (%d fine steps) exceed horizon (%d). Falling back to fine recursion.', n_fine_near, n_fine)
                end
                v_k1 = (eye(12) - F) \ C'*Qy*r_(:, k+n_preview);
                for j = n_fine:-1:1
                    v_k1 = F*v_k1 + C'*Qy*r_(:, k+j);
                end
            else
                % Decimated preview recursion
                n_fine_far   = n_fine - n_fine_near;
                n_coarse_far = floor(n_fine_far / df_(end));
                n_fine_rem   = n_fine_far - n_coarse_far * df_(end);

                tier_start = zeros(1, n_tiers);
                cursor = n_fine_rem + 1;
                for i = 1:n_tiers-1
                    tier_start(i) = cursor;
                    cursor = cursor + n_coarse_per_tier * df_(i);
                end
                tier_start(n_tiers) = cursor;  % far-horizon tier starts here

                % Precompute F^d and G_d = (I - F^d)(I - F)^{-1} = sum_{i=0}^{d-1} F^i for each tier
                G_inf = (eye(12) - F) \ eye(12);
                F_d = cell(1, n_tiers);
                G_d = cell(1, n_tiers);
                for i = 1:n_tiers
                    F_d{i} = F^df_(i);
                    G_d{i} = (eye(12) - F_d{i}) * G_inf;
                end

                % --- Backward recursion ---

                v_k1 = (eye(12) - F) \ C'*Qy*r_(:, k+n_preview);

                % (a) far-horizon tier
                for b = n_coarse_far:-1:1
                    j_blk = tier_start(n_tiers) + (b-1)*df_(end);
                    v_k1 = F_d{n_tiers}*v_k1 + G_d{n_tiers}*C'*Qy*mean(r_(:, k+j_blk : k+j_blk+df_(end)-1), 2);
                end

                % (b) near-horizon tiers
                for i = (n_tiers-1):-1:1
                    for b = n_coarse_per_tier:-1:1
                        j_blk = tier_start(i) + (b-1)*df_(i);
                        v_k1 = F_d{i}*v_k1 + G_d{i}*C'*Qy*mean(r_(:, k+j_blk : k+j_blk+df_(i)-1), 2);
                    end
                end

                % (c) remainder fine steps (computed closest to k)
                for j = n_fine_rem:-1:1
                    v_k1 = F*v_k1 + C'*Qy*r_(:, k+j);
                end
            end

            u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*v_k1);

            solve_info.NearHorizonTime = n_fine_near * qp.Ts;

        else
            % Approaching terminal state.
            % Compute full recursion for K_k, K_k^v, v_{k+1}.
            % This ensures K is ramped up appropriately according to Qyf.
            % (This is then exactly equivalent to LQT MPC using frozen A_SDC(xk), B_SDC(xk))
            P = C'*Qyf*C;
            n_preview = min(n_preview, length(r_) - k);
            v_k1 = C'*Qyf*r_(:, k+n_preview);
            n_fine = n_preview - 1;

            % Partition horizon into tiers (same structure as non-terminal branch)
            n_fine_near = n_coarse_per_tier * sum(df_(1:end-1));

            if n_fine_near >= n_fine
                % Tiers alone fill the horizon — fall back to fine recursion
                for j = n_fine:-1:1
                    K     = (R + B'*P*B) \ (B'*P*A);
                    A_cl  = A - B*K;
                    P     = C'*Qy*C + K'*R*K + A_cl'*P*A_cl;
                    v_k1  = A_cl'*v_k1 + C'*Qy*r_(:, k+j);
                end
            else
                n_fine_far   = n_fine - n_fine_near;
                n_coarse_far = floor(n_fine_far / df_(end));
                n_fine_rem   = n_fine_far - n_coarse_far * df_(end);

                % Remainder occupies j = 1..remainder, then tiers start after
                tier_start = zeros(1, n_tiers);
                cursor = n_fine_rem + 1;
                for i = 1:n_tiers-1
                    tier_start(i) = cursor;
                    cursor = cursor + n_coarse_per_tier * df_(i);
                end
                tier_start(n_tiers) = cursor;

                % Precompute coarse system matrices for each tier
                % > A_d = A^d, B_d = (I + A + ... + A^{d-1}) * B — exact
                %   ZOH discretisation at d*Ts built from the fine-step (A,B).
                % > Taking d fine steps with constant input = one coarse step.
                A_d  = cell(1, n_tiers);
                B_d  = cell(1, n_tiers);
                Qy_d = cell(1, n_tiers);
                R_d  = cell(1, n_tiers);
                A_p  = eye(12);
                B_p  = zeros(12, 6);
                p    = 0;
                for i = 1:n_tiers
                    if df_(i) == 1
                        A_d{i} = A; B_d{i} = B; Qy_d{i} = Qy; R_d{i} = R;
                    else
                        % Step up the powers to the current decimation factor
                        while p < df_(i)
                            B_p = A * B_p + B;
                            A_p = A * A_p;
                            p   = p + 1;
                        end
                        A_d{i}  = A_p;
                        B_d{i}  = B_p;
                        Qy_d{i} = df_(i) * Qy;
                        R_d{i}  = df_(i) * R;
                    end
                end

                % Backward recursion

                % (a) far-horizon tier (farthest from k)
                for b = n_coarse_far:-1:1
                    j_blk  = tier_start(n_tiers) + (b-1)*df_(end);
                    K_b    = (R_d{n_tiers} + B_d{n_tiers}'*P*B_d{n_tiers}) \ (B_d{n_tiers}'*P*A_d{n_tiers});
                    A_cl_b = A_d{n_tiers} - B_d{n_tiers}*K_b;
                    P      = C'*Qy_d{n_tiers}*C + K_b'*R_d{n_tiers}*K_b + A_cl_b'*P*A_cl_b;
                    v_k1   = A_cl_b'*v_k1 + df_(end) * C'*Qy*mean(r_(:, k+j_blk : k+j_blk+df_(end)-1), 2);
                end

                % (b) near-horizon tiers (farthest to nearest)
                for i = (n_tiers-1):-1:1
                    for b = n_coarse_per_tier:-1:1
                        j_blk  = tier_start(i) + (b-1)*df_(i);
                        K_b    = (R_d{i} + B_d{i}'*P*B_d{i}) \ (B_d{i}'*P*A_d{i});
                        A_cl_b = A_d{i} - B_d{i}*K_b;
                        P      = C'*Qy_d{i}*C + K_b'*R_d{i}*K_b + A_cl_b'*P*A_cl_b;
                        v_k1   = A_cl_b'*v_k1 + df_(i) * C'*Qy*mean(r_(:, k+j_blk : k+j_blk+df_(i)-1), 2);
                    end
                end

                % (c) remainder fine steps (nearest to k — highest value)
                for j = n_fine_rem:-1:1
                    K    = (R + B'*P*B) \ (B'*P*A);
                    A_cl = A - B*K;
                    P      = C'*Qy*C + K'*R*K + A_cl'*P*A_cl;
                    v_k1   = A_cl'*v_k1 + C'*Qy*r_(:, k+j);
                end
            end

            Kk  = (R + B'*P*B) \ (B'*P*A);
            Kvk = (R + B'*P*B) \ B';
            u   = -Kk*xk + Kvk*v_k1;

            solve_info.NearHorizonTime = n_fine_near * qp.Ts;
        end

        % One-time diagnostic: print the preview grid
        if k == 1 && ~(n_fine_near >= n_fine)
            n_coarse_tier = [repmat(n_coarse_per_tier, 1, n_tiers-1), n_coarse_far];
            n_fine_tier   = n_coarse_tier .* df_;

            fprintf('  Preview grid (%d steps/tier, %d remainder):\n', n_coarse_per_tier, n_fine_rem);
            fprintf('         DF: %s\n', sprintf('%6d', df_));
            fprintf('  FineSteps: %s  (+ %d = %d)\n', sprintf('%6d', n_fine_tier), n_fine_rem, sum(n_fine_tier) + n_fine_rem);
            fprintf('       Span: %s s  (far-horizon: %.1f Hz)\n', sprintf('%6.3f', n_fine_tier * qp.Ts), 1/(df_(end)*qp.Ts));
        end

        solve_info.DecimationFactors = df_;
        solve_info.StepsPerTier      = n_coarse_per_tier;
    end
    solve_info.TimeToComputeFeedforward = toc(clock_start);
end


function [P_ss, K_ss, info] = solve_dare(A, B, Q, R, k, solver, solverOpts)
% Solve the DARE for the steady-state Riccati solution P_ss and gain K_ss.
%
% Persistent state: warm-starts iterative solvers across calls.
% Cold-solves on k == 1 regardless of the requested solver.

    persistent P_ss_ K_ss_
    persistent warnedPreviously % for loud fallback warnings
    if k == 1, warnedPreviously = false; end

    if k == 1, solver = 'cold'; end

    switch lower(solver)

    case {'cold','sda'}
        sdaOpts = struct();
        if isfield(solverOpts, 'Tolerance')
            sdaOpts.Tolerance = solverOpts.Tolerance;
        end
        args = namedargs2cell(sdaOpts);
        [P_ss_, sinfo] = dare_sda(A, B, Q, R, args{:});
        K_ss_ = (R + B'*P_ss_*B) \ (B'*P_ss_*A);

        info.TolAchieved = sinfo.TolAchieved;
        info.NumIters    = sinfo.SolverIterations;
        info.Success     = sinfo.SolveSuccess;

    case 'idare'

        [P_ss_,K_ss_,~,sinfo] = idare(A, B, Q, R);
        parse_idare_info(sinfo);

        info.TolAchieved = compute_dare_residual(A, B, Q, R, P_ss_, K_ss_);
        info.NumIters    = 0;
        info.Success     = 1;

    case 'dlqr'

        [K_ss_,P_ss_] = dlqr(A, B, Q, R);
        info.TolAchieved = compute_dare_residual(A, B, Q, R, P_ss_, K_ss_);
        info.NumIters    = 0;
        info.Success     = 1;

    case {'nk','riccati'} % the two iterative methods

        iterOpts = solverOpts;
        iterOpts.Method = solver;
        args = namedargs2cell(iterOpts);
        [P_ss_, sinfo] = iterative_dare(A, B, Q, R, P_ss_, args{:});
        K_ss_ = (R + B'*P_ss_*B) \ (B'*P_ss_*A);

        % Fallback to SDA if initialised outside the stability basin
        % - only needed for NK
        info.Success = sinfo.SolveSuccess;
        if strcmpi(solver,'nk') && ~sinfo.SolveSuccess
            if isfield(sinfo,'UnstableK0') && sinfo.UnstableK0

                sdaOpts = struct();
                if isfield(solverOpts, 'Tolerance'), sdaOpts.Tolerance = solverOpts.Tolerance; end
                args = namedargs2cell(sdaOpts);
                [P_ss_, sinfo] = dare_sda(A, B, Q, R, args{:});
                K_ss_ = (R + B'*P_ss_*B) \ (B'*P_ss_*A);
                
                info.Success = 0.5*sinfo.SolveSuccess; % 0.5 == sentinel value for partial success
                
                if warnedPreviously == false
                    warning("solve_dare: NK iteration initialised outside of stability basin at k=%d. Fell back to SDA cold solve: success flag was %.1f", k, info.Success)
                    warnedPreviously = true;
                end
                
            end
        end

        info.TolAchieved = sinfo.TolAchieved;
        info.NumIters    = sinfo.SolverIterations;

    otherwise

        error("Unknown DARESolver type: '%s'", solver);

    end

    P_ss = P_ss_;
    K_ss = K_ss_;

end

function parse_idare_info(info)
    switch info.Report
    case 1, warning("idare(), info.Report == 1 (The solution accuracy is poor)")
    case 2, warning("idare(), info.Report == 2 (The solution is not finite)")
    case 3, error("idare(), info.Report == 3 (No solution found since the Symplectic spectrum, denoted by [L;1./L], has eigenvalues on the unit circle)")
    case 4, error("idare(), info.Report == 4 (Pencil is singular ([B;S;R] is rank deficient)")
    end
end