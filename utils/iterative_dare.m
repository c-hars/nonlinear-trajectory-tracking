function [P_ss,info] = iterative_dare(A, B, Q, R, P0, opts)
% Solve the Discrete Algebraic Riccati Equation iteratively.
%
%   P = A'PA - A'PB(R + B'PB)^{-1}B'PA + Q
%
%   Methods:
%     'riccati' - Direct DARE recursion (linear convergence, global stability from any PSD P0).
%     'nk'      - Newton-Kleinman via dlyap (default - quadratic convergence, though P0 must be stabilising).

    arguments
        A, B, Q, R, P0;
        opts.Method = 'nk';
        opts.MinItersNK = 1;
        opts.MaxItersNK = 10;  % w/ NK this can be quite low (convergence, typically within 2-4 iterations)
        opts.EarlyBreakEnabled = true;
        opts.Tolerance = 1e-4;  % DARE residual threshold. Machine precision ~= 1e-11. Can be generally be set quite a bit higher w/o perf loss. Note that - perhaps stating the obvious - tol will affect iterations/timing properties when early break is enabled. Less obviously, the tolerance in particular has a large effect on the timing variation.
        opts.RiccatiConvergenceChecksEvery = 25  % iters. Set it at 1 if you want precise breaking - when the DARE tolerance is achieved, it will break the loop and return. Set it to 25 if you're doing timing comparisons between NK and Riccati: 1 dlyap call ~= 25 Riccati recursions, in terms of compute time, and the overhead of checking the DARE tolerance (relative to one Riccati recursion) is significant, inflating the Riccati timing -- not a fair comparison if you leave it at 1. One of the nice things about Riccati is you can get fine control over the DARE tolerance, mu_dare_tol ~= p95_dare_tol, since the number of iters is so large, while NK converges in only a 3-4 iters so you have way less options re when to break the loop, can't get as consistent/repetable tolerance (just guaranteed less than spec, not that median ~= p95 like you can with Riccati). But this is just one (fairly minor) thing mind you.
    end

    if opts.MinItersNK > opts.MaxItersNK
        warning("MinItersNK (%d) exceeds MaxItersNK (%d): setting MaxItersNK = MinItersNK.", ...
                opts.MinItersNK, opts.MaxItersNK)
        opts.MaxItersNK = opts.MinItersNK;
    end

    if opts.MaxItersNK == 0
        P_ss = P0;
        info.SolverIterations = 0;
        info.TolAchieved = nan;
        info.SolveSuccess = nan;
        return
    end

    if opts.EarlyBreakEnabled && opts.MinItersNK == 0
        dare_res = compute_dare_residual(A,B,Q,R,P0);
        K    = (R + B'*P0*B) \ (B'*P0*A);
        A_cl = A - B*K;
        if (dare_res < opts.Tolerance) && all(abs(eig(A_cl)) < 1)
            P_ss = P0;
            info.SolverIterations = 0;
            info.TolAchieved = dare_res;
            info.SolveSuccess = true;
            return
        end
    end

    P = P0;

    switch lower(opts.Method)

        case 'riccati'
        % Direct iteration of the DARE map:
        %   P_{k+1} = A'P_k A - A'P_k B (R + B'P_k B)^{-1} B'P_k A + Q

            MinItersRiccati = opts.MinItersNK*25;
            MaxItersRiccati = opts.MaxItersNK*25;
            for i = 1:MaxItersRiccati
                K = (R + B'*P*B) \ (B'*P*A);
                P = A'*P*A - A'*P*B*K + Q;
                P = (P + P') / 2;

                if opts.EarlyBreakEnabled && (i >= MinItersRiccati) && ...
                        mod(i - MinItersRiccati, opts.RiccatiConvergenceChecksEvery) == 0
                    if compute_dare_residual(A,B,Q,R,P) < opts.Tolerance
                        break
                    end
                end
            end

        case 'nk'
        % Newton-Kleinman iterations solve:
        %   P_{k+1} = A_K' P_{k+1} A_K + Q + K'RK
        % via dlyap. Requires a stabilising P0. See https://arxiv.org/pdf/2503.01587.
            for i = 1:opts.MaxItersNK
                K  = (R + B'*P*B) \ (B'*P*A);
                AK = A - B*K;
                P  = dlyap(AK', Q + K'*R*K);
                P  = (P + P') / 2;

                if opts.EarlyBreakEnabled && (i >= opts.MinItersNK)
                    if compute_dare_residual(A,B,Q,R,P) < opts.Tolerance
                        break
                    end
                end
            end

        otherwise
            error('iterative_dare: Unknown method "%s". Use "riccati" or "nk".', opts.Method);
    end

    P_ss = P;
    info.SolverIterations = i;
    info.TolAchieved = compute_dare_residual(A,B,Q,R,P);
    info.SolveSuccess = (info.TolAchieved < opts.Tolerance);

end