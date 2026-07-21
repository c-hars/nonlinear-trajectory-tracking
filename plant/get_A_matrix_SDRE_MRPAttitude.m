function A = get_A_matrix_SDRE_MRPAttitude(x, qp)

    p  = x(7:9);
    p1 = p(1); p2 = p(2); p3 = p(3);
    g  = 9.81;
    A  = zeros(12);

    p_sq = p' * p;
    Sp   = S(p);
    den  = (1 + p_sq)^2;

    % Position kinematics: r_I_dot = C_Ib * v_b
    A(1:3,4:6) = eye(3) + (4*(1 - p_sq)/den)*Sp + (8/den)*(Sp*Sp);

    % Gravity SDC factorisation
    A(4:6,7:9) = 4*g / den * ...
        [ -p3, (1 - p_sq), -p1;
          (p_sq - 1),  -p3, -p2;
           2*p1,      2*p2,   0 ];

    % MRP kinematics: p_dot = B(p) * omega
    A(7:9,10:12) = 0.5 * ((1 - p_sq)/2 * eye(3) + Sp + p*p');

    % Gyroscopic coupling (Euler's rigid-body equations, SDC form)
    A(10,11) = -x(12) * (qp.I_zz - qp.I_yy) / qp.I_xx;
    A(11,12) = -x(10) * (qp.I_xx - qp.I_zz) / qp.I_yy;
    A(12,10) = -x(11) * (qp.I_yy - qp.I_xx) / qp.I_zz;
end