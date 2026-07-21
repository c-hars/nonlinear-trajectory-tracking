function [K_, xref_, ff_, Kv_, v_] = solve_LQT(r_, N, A, B, C, Qy, Qyf, R, opts)
    arguments
        r_  (:,:)   double
        N   (1,1)   {mustBeInteger}
        A   (:,:)   double
        B   (:,:)   double
        C   (:,:)   double
        Qy  (:,:)   double
        Qyf (:,:)   double
        R   (:,:)   double
        opts.UseJacobiansForRecursion = false
        opts.UseSDREFormForRecursion = false
        opts.jacobian_A_fn   function_handle
        opts.jacobian_B_fn   function_handle
        opts.SDC_A_fn   function_handle
        opts.SDC_B_fn   function_handle
        opts.qp              struct
        opts.xtraj (:,:) double
    end

    nx = size(A,1);
    nu = size(B,2);

    S_    = zeros(nx, nx, N);
    K_    = zeros(nu, nx, N);
    Kv_   = zeros(nu, nx, N);
    v_    = zeros(nx, 1, N);
    ff_   = zeros(nu, N);
    xref_ = zeros(nx, N);

    Q_lqt = C' * Qy * C;
    S_(:,:,end) = C' * Qyf * C;
    v_(:,:,end) = C' * Qyf * r_(:,N);

    for k = N-1:-1:1
        Quu = B' * S_(:,:,k+1) * B + R;
        K_(:,:,k) = -Quu \ (B' * S_(:,:,k+1) * A);
        Kk = K_(:,:,k);

        S_(:,:,k) = A' * S_(:,:,k+1) * A - Kk' * Quu * Kk + Q_lqt;
        S_(:,:,k) = 0.5 * (S_(:,:,k) + S_(:,:,k)');

        Kv_(:,:,k) = Quu \ B';
        v_(:,:,k)  = (A + B * Kk)' * v_(:,:,k+1) + C' * Qy * r_(:,k);

        % Reconstruct the state reference from the costate; S_ can be rank-deficient when the output map does not excite all states, hence the pseudoinverse
        xref_(:,k) = pinv(S_(:,:,k)) * v_(:,k);

        % Residual feedforward term
        ff_(:,k) = K_(:,:,k) * xref_(:,k) + Kv_(:,:,k) * v_(:,:,k+1);

        if opts.UseJacobiansForRecursion
            xk = opts.xtraj(:,k);
            uk = K_(:,:,k) * xk + Kv_(:,:,k) * v_(:,:,k+1);
            Ac = opts.jacobian_A_fn(xk, uk);
            Bc = opts.jacobian_B_fn(xk, uk);
            [A,B] = c2d_zoh_expm(Ac,Bc,opts.qp.Ts);
        elseif opts.UseSDREFormForRecursion
            xk = opts.xtraj(:,k);
            uk = K_(:,:,k) * xk + Kv_(:,:,k) * v_(:,:,k+1);
            Ac = opts.SDC_A_fn(xk, opts.qp);
            Bc = opts.SDC_B_fn(uk, opts.qp);
            [A,B] = c2d_zoh_expm(Ac,Bc,opts.qp.Ts);
        end
    end

    xref_(:,end) = pinv(S_(:,:,end)) * v_(:,end);
end