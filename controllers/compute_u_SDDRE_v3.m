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
        opts.DecimationFactor = 1 % scalar or ascending vector of DFs (geometric tiers). Tiers 1..end-1 are intermediate (each gets StepsPerTier decimated steps); last element is the tail DF filling the rest of the horizon.
        opts.StepsPerTier     = 3 % number of decimated steps per intermediate tier
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

    clock_start = tic;
    Q = C'*Qy*C;
    [P_ss, K_ss, dare_info] = solve_dare(A, B, Q, R, k, opts.DARESolver, opts.DARESolverOpts);
    solve_info.DARESolverTolAchieved = dare_info.TolAchieved;
    solve_info.DARESolverNumIters    = dare_info.NumIters;
    solve_info.DARESolverSuccess     = dare_info.Success;
    solve_info.TimeToSolveDARE = toc(clock_start);
    

    % --- 3. Compute the feedforward (reference preview) term ---
    
    clock_start = tic;
    A_cl  = A - B*K_ss;

    persistent CtQyr_
    % Weighted reference, C'*Qy*r_
    % This quantity is referenced throughout - for efficiency, prefer precomputing as a fixed offline table. Recomputing it each call is redundant computation.
    % But if the reference plan changes, or you're feeding a constantly changing reference, not a known-apriori plan - then you can't get away with precomputing. A circular buffer sized to the preview window (and only updated with the one-step new additions) would be a good idea there.
    % Padded by M+d_max so no index clamping is needed
    if k == 1 || isempty(CtQyr_)
        d_max = max(opts.DecimationFactor);
        M = round(opts.PreviewHorizon / qp.Ts);
        rpad = [r_, repmat(r_(:,end), 1, M + d_max)];
        CtQyr_ = C' * Qy * rpad;
    end

    if isinf(opts.PreviewHorizon)
        % Constant reference approximation
        s_k1  = (eye(12) - A_cl') \ (C'*Qy*r_(:,k));
        u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*s_k1);

    else
        M = round(opts.PreviewHorizon/qp.Ts);   % receding horizon preview window
        approachingTerminal = (k+M >= length(r_));

        % Geometric tier decimation setup (shared by both branches)
        df = sort(opts.DecimationFactor(:)', 'ascend');
        nb = opts.StepsPerTier;
        ntiers = length(df);

        if ~(opts.AlwaysUseFullFiniteHorizonMPC || (opts.UseFullFiniteHorizonMPCAtTerminal && approachingTerminal))

            F = A_cl';
            nsteps = M - 1;

            % Partition horizon into geometric tiers
            %   Layout from k outward (j=1 nearest):
            %     leftover | tier 1 (nb blocks @ df(1)) | tier 2 (nb @ df(2)) | ... | tail (df(end))
            %   Leftover fine steps sit at the near end — highest marginal value for resolution.
            nearfield_fine = nb * sum(df(1:end-1));  % fine steps consumed by intermediate tiers

            if nearfield_fine >= nsteps
                % Tiers alone fill the horizon — fall back to full-rate recursion
                if k == 1
                    warning('compute_u_SDDRE_v3: nearfield tiers (%d fine steps) exceed horizon (%d). Falling back to full-rate.', nearfield_fine, nsteps)
                end
                v_k1 = (eye(12) - F) \ CtQyr_(:, k+M);
                for j = nsteps:-1:1
                    v_k1 = F*v_k1 + CtQyr_(:, k+j);
                end
            else
                d_tail = df(end);
                tail_avail = nsteps - nearfield_fine;
                tail_nblk  = floor(tail_avail / d_tail);
                leftover   = nsteps - nearfield_fine - tail_nblk * d_tail;

                % Tier start indices (fine-step offset from k+1, nearest first)
                % Leftover occupies j = 1..leftover, then tiers start after
                tier_lo = zeros(1, ntiers);
                cursor = leftover + 1;
                for i = 1:ntiers-1
                    tier_lo(i) = cursor;
                    cursor = cursor + nb * df(i);
                end
                tier_lo(ntiers) = cursor;  % tail starts here

                % Precompute F^d and G_d = (I - F^d)(I - F)^{-1} = sum_{i=0}^{d-1} F^i for each tier
                Ginf = (eye(12) - F) \ eye(12);
                Fd = cell(1, ntiers);
                Gd = cell(1, ntiers);
                for i = 1:ntiers
                    Fd{i} = F^df(i);
                    Gd{i} = (eye(12) - Fd{i}) * Ginf;
                end

                % Backward recursion
                v_k1 = (eye(12) - F) \ CtQyr_(:, k+M);

                % (a) tail blocks (farthest from k)
                for b = tail_nblk:-1:1
                    jblk = tier_lo(ntiers) + (b-1)*d_tail;
                    v_k1 = Fd{ntiers}*v_k1 + Gd{ntiers}*mean(CtQyr_(:, k+jblk : k+jblk+d_tail-1), 2);
                end

                % (b) intermediate tiers (farthest to nearest)
                for i = (ntiers-1):-1:1
                    d_i = df(i);
                    for b = nb:-1:1
                        jblk = tier_lo(i) + (b-1)*d_i;
                        v_k1 = Fd{i}*v_k1 + Gd{i}*mean(CtQyr_(:, k+jblk : k+jblk+d_i-1), 2);
                    end
                end

                % (c) leftover fine steps (nearest to k — highest value)
                for j = leftover:-1:1
                    v_k1 = F*v_k1 + CtQyr_(:, k+j);
                end
            end

            u = -K_ss*xk + (R + B'*P_ss*B) \ (B'*v_k1);

            solve_info.NearfieldTime = nearfield_fine * qp.Ts;

        else
            % Approaching terminal state.
            % Compute full recursion for K_k, K_k^v, v_{k+1}.
            % This ensures K is ramped up appropriately according to Qyf.
            % (This is then exactly equivalent to LQT MPC using frozen A_SDC(xk), B_SDC(xk))
            P = C'*Qyf*C;
            M = min(M, length(r_) - k);
            v = C'*Qyf*r_(:, k+M);
            nsteps = M - 1;

            % Partition horizon into geometric tiers (same structure as non-terminal branch)
            nearfield_fine = nb * sum(df(1:end-1));

            if nearfield_fine >= nsteps
                % Tiers alone fill the horizon — fall back to full-rate recursion
                for j = nsteps:-1:1
                    K_j    = (R + B'*P*B) \ (B'*P*A);
                    A_cl_j = A - B*K_j;
                    P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                    v      = A_cl_j'*v + CtQyr_(:, k+j);
                end
            else
                d_tail = df(end);
                tail_avail = nsteps - nearfield_fine;
                tail_nblk  = floor(tail_avail / d_tail);
                leftover   = nsteps - nearfield_fine - tail_nblk * d_tail;

                % Leftover occupies j = 1..leftover, then tiers start after
                tier_lo = zeros(1, ntiers);
                cursor = leftover + 1;
                for i = 1:ntiers-1
                    tier_lo(i) = cursor;
                    cursor = cursor + nb * df(i);
                end
                tier_lo(ntiers) = cursor;

                % Precompute decimated system matrices for each tier
                % > Exact ZOH discretisation at d*Ts, built from the fine-step (A,B)  
                % > Taking d fine steps with a constant input is the same as one coarse step
                %   with Ad = A^d and Bd = (I + A + A^2 + ... + A^{d-1}) * B = Gd * B
                Adt  = cell(1, ntiers);
                Bdt  = cell(1, ntiers);
                Qydt = cell(1, ntiers);
                Rdt  = cell(1, ntiers);
                for i = 1:ntiers
                    d = df(i);
                    if d == 1
                        Adt{i} = A; Bdt{i} = B; Qydt{i} = Qy; Rdt{i} = R;
                    else
                        Apow = eye(12);
                        G    = eye(12);
                        for p = 1:d-1
                            Apow = A * Apow;
                            G    = G + Apow;
                        end
                        Adt{i}  = A * Apow;
                        Bdt{i}  = G * B;
                        Qydt{i} = Qy * d;
                        Rdt{i}  = R * d;
                    end
                end

                % Backward recursion

                % (a) tail blocks (farthest from k)
                for b = tail_nblk:-1:1
                    jblk   = tier_lo(ntiers) + (b-1)*d_tail;
                    K_b    = (Rdt{ntiers} + Bdt{ntiers}'*P*Bdt{ntiers}) \ (Bdt{ntiers}'*P*Adt{ntiers});
                    A_cl_b = Adt{ntiers} - Bdt{ntiers}*K_b;
                    P      = C'*Qydt{ntiers}*C + K_b'*Rdt{ntiers}*K_b + A_cl_b'*P*A_cl_b;
                    v      = A_cl_b'*v + d_tail * mean(CtQyr_(:, k+jblk : k+jblk+d_tail-1), 2);
                end

                % (b) intermediate tiers (farthest to nearest)
                for i = (ntiers-1):-1:1
                    d_i = df(i);
                    for b = nb:-1:1
                        jblk   = tier_lo(i) + (b-1)*d_i;
                        K_b    = (Rdt{i} + Bdt{i}'*P*Bdt{i}) \ (Bdt{i}'*P*Adt{i});
                        A_cl_b = Adt{i} - Bdt{i}*K_b;
                        P      = C'*Qydt{i}*C + K_b'*Rdt{i}*K_b + A_cl_b'*P*A_cl_b;
                        v      = A_cl_b'*v + d_i * mean(CtQyr_(:, k+jblk : k+jblk+d_i-1), 2);
                    end
                end

                % (c) leftover fine steps (nearest to k — highest value)
                for j = leftover:-1:1
                    K_j    = (R + B'*P*B) \ (B'*P*A);
                    A_cl_j = A - B*K_j;
                    P      = C'*Qy*C + K_j'*R*K_j + A_cl_j'*P*A_cl_j;
                    v      = A_cl_j'*v + CtQyr_(:, k+j);
                end
            end

            Kk  = (R + B'*P*B) \ (B'*P*A);
            Kvk = (R + B'*P*B) \ B';
            u   = -Kk*xk + Kvk*v;

            solve_info.NearfieldTime = nearfield_fine * qp.Ts;
        end

        % One-time diagnostic: print the decimation schedule
        if k == 1 && ~(nearfield_fine >= nsteps)
            tail_nblk = floor((M - 1 - nearfield_fine) / df(end));
            nblks = [repmat(nb, 1, ntiers-1), tail_nblk];
            fine  = nblks .* df;

            fprintf('  Decimation schedule (%d steps/tier, %d leftover):\n', nb, leftover);
            fprintf('         DF: %s\n', sprintf('%6d', df));
            fprintf('  FineSteps: %s  (+ %d = %d)\n', sprintf('%6d', fine), leftover, sum(fine) + leftover);
            fprintf('       Span: %s s  (tail: %.1f Hz)\n', sprintf('%6.3f', fine * qp.Ts), 1/(df(end)*qp.Ts));
        end

        solve_info.DecimationTiers = df;
        solve_info.StepsPerTier    = nb;
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