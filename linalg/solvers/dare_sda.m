function [P, info] = dare_sda(A, B, Q, R, opts)
% Solves the discrete-time algebraic Riccati equation:
%   P = A'PA - A'PB(R + B'PB)^{-1}B'PA + Q
% 
% Same calling convention as MATLAB-native dare(); drop-in replacement.
% 
% Fast, portable cold solver that can be implemented on embedded. Faster *suboptimal solves* (via reduced tol) are also possible. 
% 
% Uses the structure-preserving doubling algorithm (SDA):
%       Wk      = I + G_k*H_k
%       A_{k+1} = A_k*(Wk\A_k)
%       G_{k+1} = G_k + A_k*(Wk\G_k)*A_k'
%       H_{k+1} = H_k + A_k'*H_k*(Wk\A_k)
% (doi:10.1080/00207170410001714988, doi:10.1002/gamm.202000018)
%
% Preconditions are unchecked by design (this is a hot-loop solver).

    arguments
        A, B, Q, R;
        opts.Tolerance    = 1e-4;   % DARE residual threshold
        opts.MinDoublings = 1;
        opts.MaxDoublings = 40;     % covers any reasonable rho (up to ~1-1e-9)
    end

    n  = size(A,1);
    In = eye(n);

    Ak = A;
    G  = B*(R\B');  G = (G + G')/2;
    H  = Q;

    for k = 1:opts.MaxDoublings
        Wk = In + G*H;

        S  = Wk \ [Ak, G];  % one factorisation, two right-hand sides
        Yk = S(:, 1:n);
        Zk = S(:, n+1:end);

        dH = Ak'*H*Yk;   H  = H + dH;  H = (H + H')/2;
        dG = Ak*Zk*Ak';  G  = G + dG;  G = (G + G')/2;
        Ak = Ak*Yk;

        if k >= opts.MinDoublings
            % norm(Ak) is the algorithm's native monitor and nearly free. It also catches the case where further doublings cannot progress, which the residual test alone would not.
            if norm(Ak, 'fro') < eps
                break
            end
            if compute_dare_residual(A,B,Q,R,H) < opts.Tolerance
                break  % norm(dH) is the cheaper check (used in the original publication)
            end
        end
    end

    P = H;
    info.SolverIterations = k;
    info.TolAchieved      = compute_dare_residual(A,B,Q,R,P);
    info.SolveSuccess     = info.TolAchieved < opts.Tolerance;
end 
