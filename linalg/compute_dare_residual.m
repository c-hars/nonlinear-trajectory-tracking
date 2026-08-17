function dare_res = compute_dare_residual(A,B,Q,R,P,K)
    if nargin < 6
        K = (R + B'*P*B) \ (B'*P*A);
    end
    A_cl = A - B*K;
    dare_res = norm(A_cl'*P*A_cl - P + Q + K'*R*K, 'fro') / norm(P);
    % den1 = norm(P, 'fro');
    % den2 = norm(P, 'fro') + norm(Q, 'fro') + norm(K'*(R+B'*P*B)*K, 'fro');
    % fprintf("den1 = %.3e, den2 = %.3e, err=%.3e\n", den1,den2,abs(den2-den1));
end