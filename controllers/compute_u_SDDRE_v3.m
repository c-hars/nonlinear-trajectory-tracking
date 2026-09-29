function [u,solve_info] = compute_u_SDDRE_v3(tk,xk,k,uk,r_,C,Qy,R,Qyf,qp,opts)
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
        qp
        opts.SDC_A_function = @ get_A_matrix_SDRE_EulerAttitude
        opts.SDC_B_function = @ get_B_matrix_SDRE
        opts.PreviewHorizon = 2.0 % this is all really designed for finite-horizon formulations, with the horizon large (>= 2.0 for our case) -- but inf is a safe choice, just get delayed tracking
        opts.UseFullFiniteHorizonMPCAtTerminal = true % best set to false if you need precisely consistent/predictable solve times - adds a bit of overhead to do the full recursion for K too
        opts.AlwaysUseFullFiniteHorizonMPC = false % set to true → regression to the standard MPC cost function being optimised / directly comparable OCP with LinMPC, NLMPC -- no longer optimal preview control that uses the infinite/finite horizon split (core part of the algorithm!)
        opts.DecimationOpts = struct('DecimationFactors', [1], 'StepsPerTier', 1)

        opts.DARESolver (1,:) char {mustBeMember(opts.DARESolver, {'nk','riccati','cold','idare','dlqr','sda'})} = 'nk'
        opts.DARESolverOpts (1,1) struct = struct()
    end

    % --- 0. Bookkeeping ---
    % run_sim() needs homogenous solve_info: populate upfront with dummy data
    solve_info = struct( ...
        'TimeToComputeDiscreteSDCMatrices', 0, ...
        'DARESolverTolAchieved', NaN, ...
        'DARESolverNumIters', NaN, ...
        'DARESolverSuccess', NaN, ...
        'TimeToSolveDARE', 0, ...
        'TimeToComputeFeedforward', 0);


    % --- 1. Compute the (ZOH-discretised) SDC matrices ---

    clock_start = tic;
    Ac = opts.SDC_A_function(xk, qp);
    uk = max(uk, -qp.nominal_omegas(:));
    uk = min(uk,  qp.max_du(:));
    Bc = opts.SDC_B_function(uk, qp);
    [A,B] = c2d_zoh_expm(Ac,Bc,qp.Ts);
    solve_info.TimeToComputeDiscreteSDCMatrices = toc(clock_start);


    % --- 2. Compute the control inputs: u = -K*x + Kv*v ---

    n_preview = round(opts.PreviewHorizon/qp.Ts);

    % One-time diagnostic: print the preview grid
    if k == 1 && ~isinf(n_preview)
        print_horizon_grid(build_horizon_grid(n_preview, opts.DecimationOpts), qp);
    end

    approachingTerminal = ~isinf(n_preview) && (k + n_preview >= length(r_));
    using_full_finite_horizon_mpc = opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal);
    using_preview_control = ~using_full_finite_horizon_mpc;
    
    if using_preview_control

    % --- 2.1: SD-OPT path ---
    % Preview control: uses the infinite horizon gain/cost-to-go plus a reference preview term.

        % --- 2.1.1: Solve the DARE for P_ss, K_ss ---

        clock_start = tic;
        Q = C'*Qy*C;
        [P_ss, K_ss, dare_info] = solve_dare(A, B, Q, R, k, opts.DARESolver, opts.DARESolverOpts);
        solve_info.DARESolverTolAchieved = dare_info.TolAchieved;
        solve_info.DARESolverNumIters    = dare_info.NumIters;
        solve_info.DARESolverSuccess     = dare_info.Success;
        solve_info.TimeToSolveDARE = toc(clock_start);

        % --- 2.1.2: Compute the preview term v_{k+1} ---

        clock_start = tic;
        A_cl  = A - B*K_ss;
        if isinf(opts.PreviewHorizon)
            % Constant reference approximation: no recursion
            v_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        else
            % Recursion over the preview window
            r_ = r_(:, min(1:k+n_preview, end));  % pad the reference if needed
            F = A_cl';
            if isempty(opts.DecimationOpts) || all(opts.DecimationOpts.DecimationFactors == 1)
                % No decimation: standard Riccati recursion for the preview
                v_k1 = (eye(12) - F) \ C'*Qy*r_(:, k+n_preview);  % SDOPT seed
                for j = n_preview-1:-1:1
                    v_k1 = F*v_k1 + C'*Qy*r_(:, k+j);
                end
            else
                % Decimation: the preview grid is constructed and processed on a coarser layout - not all at the full sample rate - as determined by DecimationOpts
                grid = build_horizon_grid(n_preview, opts.DecimationOpts);
                v_k1 = compute_preview(F, C, Qy, r_, k, n_preview, grid);
            end
        end
        solve_info.TimeToComputeFeedforward = toc(clock_start);

        % --- 2.1.3: Control input ---

        Kv = (R + B'*P_ss*B) \ B';
        u = -K_ss*xk + Kv*v_k1;


    else

    % --- 2.2: SD-MPC path ---
    % State-dependent MPC: runs the full Riccati recursion to get both the gain and preview terms. No DARE
    % This branch is primarily used for shrinking horizon MPC (at the terminal stage)
    % Also can be used to run standard MPC machinery on the whole problem, not just the terminal part.

        clock_start = tic;

        % --- 2.2.1: Horizon and terminal conditions ---

        n_preview = min(n_preview, length(r_) - k);  % shrinks as the end of the reference approaches
        P    = C'*Qyf*C;
        v_k1 = C'*Qyf*r_(:, k+n_preview);

        % --- 2.2.2: Riccati recursion (gain + preview) ---

        if isempty(opts.DecimationOpts) || all(opts.DecimationOpts.DecimationFactors == 1)
            % No decimation: standard Riccati recursion
            for j = n_preview-1:-1:1
                K    = (R + B'*P*B) \ (B'*P*A);
                A_cl = A - B*K;
                P    = C'*Qy*C + K'*R*K + A_cl'*P*A_cl;
                v_k1 = A_cl'*v_k1 + C'*Qy*r_(:, k+j);
            end
        else
            grid = build_horizon_grid(n_preview, opts.DecimationOpts);
            [P, v_k1] = compute_riccati_recursion(A, B, C, Qy, R, P, v_k1, r_, k, grid);
        end

        % --- 2.2.3: Control input ---

        solve_info.TimeToComputeFeedforward = toc(clock_start); % Note: slight misnomer here - should be TimeToComputeRecursion. Kept for now (maintains the same struct format as the SDOPT branch).

        Kk  = (R + B'*P*B) \ (B'*P*A);
        Kvk = (R + B'*P*B) \ B';
        u   = -Kk*xk + Kvk*v_k1;

    end

end

function print_horizon_grid(grid, qp)
% One-time diagnostic: print the preview grid

    if grid.use_fine_fallback, return, end

    n_fine_tier = grid.n_blocks .* grid.df_;

    fprintf('  Preview grid (%d steps/tier, %d remainder):\n', grid.n_steps_per_tier, grid.n_fine_rem);
    fprintf('         DF: %s\n', sprintf('%6d', grid.df_));
    fprintf('  FineSteps: %s  (+ %d = %d)\n', sprintf('%6d', n_fine_tier), grid.n_fine_rem, sum(n_fine_tier) + grid.n_fine_rem);
    fprintf('       Span: %s s  (far-horizon: %.1f Hz)\n', sprintf('%6.3f', n_fine_tier * qp.Ts), 1/(grid.df_(end)*qp.Ts));
end

function grid = build_horizon_grid(n_preview, DecimationOpts)
% Partition the preview horizon into remainder fine steps, near-horizon tiers and the far-horizon tier.

    df_ = sort(DecimationOpts.DecimationFactors(:)', 'ascend');
    n_steps_per_tier = DecimationOpts.StepsPerTier;
    n_tiers = length(df_);

    n_fine      = n_preview - 1;
    n_fine_near = n_steps_per_tier * sum(df_(1:end-1));

    grid.df_               = df_;
    grid.n_tiers           = n_tiers;
    grid.n_steps_per_tier  = n_steps_per_tier;
    grid.n_fine            = n_fine;
    grid.n_fine_near       = n_fine_near;
    grid.use_fine_fallback = (n_fine_near >= n_fine);  % tiers alone fill the horizon

    if grid.use_fine_fallback
        return
    end

    n_fine_far   = n_fine - n_fine_near;
    n_coarse_far = floor(n_fine_far / df_(end));
    n_fine_rem   = n_fine_far - n_coarse_far * df_(end);
    grid.n_blocks = [repmat(n_steps_per_tier, 1, n_tiers-1), n_coarse_far];  % coarse steps per tier

    tier_start = zeros(1, n_tiers);
    cursor = n_fine_rem + 1;
    for i = 1:n_tiers-1
        tier_start(i) = cursor;
        cursor = cursor + n_steps_per_tier * df_(i);
    end
    tier_start(n_tiers) = cursor;  % far-horizon tier starts here

    grid.n_coarse_far = n_coarse_far;
    grid.n_fine_rem   = n_fine_rem;
    grid.tier_start   = tier_start;
end

function v_k1 = compute_preview(F, C, Qy, r_, k, n_preview, grid)
% Decimated Riccati recursion (costate only - under preview control the cost-to-go and gain are fixed via P=P_ss, K=K_ss)

    v_k1 = (eye(12) - F) \ C'*Qy*r_(:, k+n_preview);  % SDOPT seed

    if grid.use_fine_fallback
        % Fall back to fine recursion
        for j = grid.n_fine:-1:1
            v_k1 = F*v_k1 + C'*Qy*r_(:, k+j);
        end
        return
    end

    df_ = grid.df_;
    n_tiers = grid.n_tiers;

    % Precompute the coarse preview matrices for each tier.
    % > The fine preview update is:
    %       v <- F*v + C'*Qy*r
    %   With r held constant, d fine steps collapse into one coarse step:
    %       v <- F_d*v + G_d*C'*Qy*r
    %   where
    %       F_d = F^d
    %   and
    %       G_d = I + F + ... + F^{d-1}.
    % > G_d is a geometric series, which can be efficiently evaluated via
    %       G_d = (I - F^d)(I - F)^{-1}
    %   since F^d is already computed.
    G_inf = (eye(12) - F) \ eye(12);
    F_d = cell(1, n_tiers);
    G_d = cell(1, n_tiers);
    for i = 1:n_tiers
        F_d{i} = F^df_(i);
        G_d{i} = (eye(12) - F_d{i}) * G_inf;
    end

    % Recursion across the tiers (far-horizon and near-horizon)
    for i = n_tiers:-1:1
        d = df_(i);
        for b = grid.n_blocks(i):-1:1
            j_blk = grid.tier_start(i) + (b-1)*d;
            r_bar = sum(r_(:, k+j_blk : k+j_blk+d-1), 2) / d;  % block-averaged reference
            v_k1  = F_d{i}*v_k1 + G_d{i}*C'*Qy*r_bar;
        end
    end

    % Recursion across the remainder fine steps
    for j = grid.n_fine_rem:-1:1
        v_k1 = F*v_k1 + C'*Qy*r_(:, k+j);
    end
end

function [P, v_k1] = compute_riccati_recursion(A, B, C, Qy, R, P, v_k1, r_, k, grid)
% Decimated Riccati recursion (cost-to-go and costate)

    Q = C'*Qy*C;

    if grid.use_fine_fallback
        % Tiers alone fill the horizon - fall back to fine recursion
        for j = grid.n_fine:-1:1
            K    = (R + B'*P*B) \ (B'*P*A);
            A_cl = A - B*K;
            P    = Q + K'*R*K + A_cl'*P*A_cl;
            v_k1 = A_cl'*v_k1 + C'*Qy*r_(:, k+j);
        end
        return
    end

    df_ = grid.df_;
    n_tiers = grid.n_tiers;

    % Precompute the coarse system matrices for each tier.
    % > Since the decimation factors are sorted in ascending order, the powers A^p and B_p can be built up incrementally as the decimation factor increases (rather than recomputed from scratch each time).
    A_d  = cell(1, n_tiers);
    B_d  = cell(1, n_tiers);
    Qy_d = cell(1, n_tiers);
    Q_d  = cell(1, n_tiers);
    R_d  = cell(1, n_tiers);
    A_p  = eye(12);
    B_p  = zeros(12, 6);
    p    = 0;
    for i = 1:n_tiers
        if df_(i) == 1
            A_d{i} = A; B_d{i} = B; Qy_d{i} = Qy; Q_d{i} = Q; R_d{i} = R;
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
            Q_d{i}  = df_(i) * Q;
            R_d{i}  = df_(i) * R;
        end
    end

    % Recursion across the tiers (far-horizon and near-horizon)
    for i = n_tiers:-1:1
        d = df_(i);
        for b = grid.n_blocks(i):-1:1
            j_blk = grid.tier_start(i) + (b-1)*d;
            r_bar = sum(r_(:, k+j_blk : k+j_blk+d-1), 2) / d;
            K_d   = (R_d{i} + B_d{i}'*P*B_d{i}) \ (B_d{i}'*P*A_d{i});
            A_cl  = A_d{i} - B_d{i}*K_d;
            P     = Q_d{i} + K_d'*R_d{i}*K_d + A_cl'*P*A_cl;
            v_k1  = A_cl'*v_k1 + C'*Qy_d{i}*r_bar;
        end
    end

    % Recursion across the remainder fine steps
    for j = grid.n_fine_rem:-1:1
        K    = (R + B'*P*B) \ (B'*P*A);
        A_cl = A - B*K;
        P    = Q + K'*R*K + A_cl'*P*A_cl;
        v_k1 = A_cl'*v_k1 + C'*Qy*r_(:, k+j);
    end
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