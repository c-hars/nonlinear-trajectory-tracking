function [X, info] = dlyap_sda(A, Q, tol)
% Solves the discrete-time Lyapunov equation:
%   A*X*A' - X + Q = 0
% 
% Same calling convention as MATLAB-native dlyap(); drop-in replacement.
% 
% [X, info] = dlyap_sda(A, Q) additionally returns an info struct; the key piece is info.IsStable, a Schur stability verdict on A. Used in iterative_dare (NK iteration cannot begin from an unstable initial guess).
% 
% Uses Smith squaring / a structure-preserving doubling algorithm (SDA); readily portable to embedded, and also is iterative - allows solves down to machine precision, or to some prespecified tolerance (early break), or simply for some fixed budget of iterations.
% 
% Unlike dlyap or dlyap_schur, a closed-loop stable A *is* required to solve the Lyapunov equation via this method. Failure to converge within the iteration limit yields an instability verdict for A.
% 
% Preconditions (compatible sizes, Q positive semi-definite, (A,Q^{1/2}) observable) are unchecked by design (this is a hot-loop solver).
% 
% References:
%   doi:10.1137/0116017 (introduction of the Smith doubling technique, originally applied to the continuous time Sylvester equation),
%   doi:10.1002/gamm.202000018 (modern overview of SDA algorithms, including the dlyap (Stein equation) case).

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