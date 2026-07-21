function A = get_A_matrix_SDRE_FRAAttitude(x, qp)

    th  = x(7:9);
    g   = 9.81;
    A   = zeros(12);

    Sth  = S(th);
    ang  = norm(th);
    g_   = sinc1(ang);
    f_   = 0.5 * sinc1(ang/2)^2;
    if ang < 1e-2
        h_ = (1/6) - ang^2/120 + ang^4/5040;
    else
        h_ = (ang - sin(ang)) / ang^3;
    end

    % Position kinematics: r_I_dot = C_Ib * v_b
    A(1:3,4:6) = eye(3) + g_*Sth + f_*(Sth*Sth);

    % Gravity SDC factorisation with symmetric split on cross terms
    th1 = th(1); th2 = th(2); th3 = th(3);
    A(4:6,7:9) = g * [ -f_*th3/2,        g_,  -f_*th1/2;
                             -g_,  -f_*th3/2,  -f_*th2/2;
                          f_*th1,     f_*th2,          0 ];

    % FRA kinematics: theta_dot = inv(F) * omega
    F = eye(3) - f_*Sth + h_*(Sth*Sth);
    A(7:9,10:12) = inv(F);

    % Gyroscopic coupling (Euler's rigid-body equations, SDC form)
    A(10,11) = -x(12) * (qp.I_zz - qp.I_yy) / qp.I_xx;
    A(11,12) = -x(10) * (qp.I_xx - qp.I_zz) / qp.I_yy;
    A(12,10) = -x(11) * (qp.I_yy - qp.I_xx) / qp.I_zz;
end