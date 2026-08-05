function [X, info] = dlyap_sda(A, Q, tol)
% ---
% 
% Solves the discrete-time Lyapunov equation:
%   A*X*A' - X + Q = 0
% 
% Same calling convention as MATLAB-native dlyap(); drop-in replacement.
% Also returns an info struct (notably info.IsStable: Schur stability verdict on A).
% 
% Fast, portable cold solver. Faster *suboptimal solves* (via reduced tol) are also possible.
% Preconditions are unchecked by design (hot-loop solver).
% Convergence requires a Schur-stable A (providing a stability test).
% 
% ---
% 
% Smith doubling evaluates the series:
%       X = sum_{i=0}^{inf} A^i * Q * (A')^i
% via the partial-sum doubling identity:
%       S_{2N} = S_N + A^N * S_N * (A^N)'
% where S_N = sum_{i=0}^{N-1} A^i * Q * (A')^i (i.e. the N-truncated sum).
% Iterations thus yield S_1 -> S_2 -> S_4 -> S_8 -> ... (quadratic convergence).
% 
% In the implementation: X tracks S_N, P tracks A^N. Each iteration computes:
%       dX = P*X*P', X <- X + dX, P <- P^2
% yielding X = S_2, S_4, S_8, ... for j = 1, 2, 3, ... .
% 
% ---
% 
% (doi:10.1137/0116017, doi:10.1002/gamm.202000018)

    cls = class(A);
    if nargin < 3 || isempty(tol)
        if strcmp(cls, 'single')
            tol = 1e-6;
        else
            tol = 1e-14;
        end
    end

    if strcmp(cls, 'single')
        maxDoublings = 25; % Practical single-precision cutoff; reliably covers rho up to ~1-16*eps, eps = 1.2e-07 (= 2^-23)
    else
        maxDoublings = 55; % Practical double-precision cutoff; reliably covers rho up to ~1-16*eps, eps = 2.2e-16 (= 2^-52)
    end
    epsC = eps(cls);
    info.Converged = false;
    info.IsStable  = false;
    ninc = 0;
    ninc_prev = 0;

    X = Q;
    P = A;
    for j = 1:maxDoublings
        inc  = P * X * P';
        X    = X + inc;
        ninc = norm(inc, 'fro');
        nX   = norm(X, 'fro');
        if j > 1 && ninc_prev > 0
            r = ninc / ninc_prev;
            if r < 1 && ninc * r / (1 - r) <= tol * nX
            % NB: ninc * r / (1 − r) is a (conservative) geometric extrapolation of the remaining tail
            % We've converged if the remaining error is below tolerance
                info.Converged = true;
                info.IsStable  = true;
                break
            end
        end
        ninc_prev = ninc;
        P = P * P;
    end

    if ishermitian(Q)
        X = (X + X') / 2;
    end

    info.Doublings    = j;
    if nargout > 1 % relatively expensive computations: gated by necessity
        info.RelIncrement = ninc / max(norm(X, 'fro'), epsC);
        info.Residual = norm(A*X*A' - X + Q, 'fro') / max(norm(Q, 'fro'), epsC);
    end
end