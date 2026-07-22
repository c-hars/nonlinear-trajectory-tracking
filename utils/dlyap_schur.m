function [X, info] = dlyap_schur(A, Q)
% Solves the discrete-time Lyapunov equation:
%   A*X*A' - X + Q = 0
% 
% Same calling convention as MATLAB-native dlyap(); drop-in replacement.
%
% [X, info] = dlyap_schur(A, Q) additionally returns an info struct; the key piece is info.IsStable, which is a stability check on A (Schur stability, i.e. SpectralRadius < 1).
% 
% Schur-stable A is not required to solve the Lyapunov equation (only non-uniqueness); MATLAB-native dlyap will happily return the (algebraically valid) answer for unstable closed-loop matrices. But SpectralRadius < 1 is required for a meaningful solution in our usage - starting NK iteration from an unstable closed-loop matrix A_cl = A-B*K is not valid; this function is essentially an exact copy of MATLAB's native dlyap, but with that additional stability check obtained for free (we have the eigenvalues from the Schur decomposition stage).
%
% Preconditions are unchecked by design (this is a hot-loop solver).

    n = size(A, 1);

    % Schur transform
    [U, T] = schur(A);
    Qt = U' * Q * U;

    % Solve Y - T*Y*T' = Qt in one compiled call
    % This is the same call that MATLAB's dlyap makes
    % Note that 'transp' selects the ...*T' orientation
    Y = matlab.internal.math.sylvester_tri(T, 'I', Qt, 'I', -T, 'transp');

    % Transform back and symmetrise away roundoff
    X = U * Y * U';
    X = (X + X') / 2;

    % Certificate (only paid for when asked for)
    if nargout > 1
        ev = ev_from_blocks(T, n);
        info.Eigenvalues    = ev;
        info.SpectralRadius = max(abs(ev));
        info.IsStable       = info.SpectralRadius < 1;
        info.MinDivisor     = min(abs(1 - ev*ev'), [], 'all'); % divisors are 1 - lambda_i*conj(lambda_j)
        info.Residual       = norm(A*X*A' - X + Q, 'fro') / ...
                              max(norm(Q, 'fro'), eps);
    end
end

function ev = ev_from_blocks(T, n)
% Eigenvalues from the diagonal blocks of a (quasi) triangular T.
% A nonzero subdiagonal entry marks a 2x2 block. Note that LAPACK writes exact zeros between blocks, so no tolerance is needed.
    ev = diag(T);
    if n < 2 || ~isreal(T)
        return  % complex Schur: diagonal is exact
    end
    k = find(diag(T, -1) ~= 0);  % row index of each 2x2 block start
    if isempty(k)
        return
    end
    sup = T(k*n + k);
    sub = T((k-1)*n + k+1);
    m   = (ev(k) + ev(k+1)) / 2;
    s   = sqrt(complex(m.^2 - (ev(k).*ev(k+1) - sup.*sub)));
    ev(k)   = m + s;
    ev(k+1) = m - s;
end