function dare_res = compute_dare_residual(A,B,Q,R,P,K)
    if nargin < 6
        K = (R + B'*P*B) \ (B'*P*A);
    end
    A_cl = A - B*K;
    dare_res = norm(A_cl'*P*A_cl - P + Q + K'*R*K, 'fro');
end