function Ac = get_A_matrix_SDRE_EulerAttitude(x, qp)

    phi = x(7); theta = x(8); psi = x(9);

    Ac = get_A_matrix();

    % Position kinematics: r_I_dot = C_Ib * v_b
    C_bI = C_x(phi) * C_y(theta) * C_z(psi);
    Ac(1:3,4:6) = C_bI';

    % Gravity SDC factorisation of C_bI*[0;0;-g] residual about hover
    g = 9.81;
    Ac(4,8) =  g * sinc1(theta);
    Ac(5,7) = -g * cos(theta) * sinc1(phi);
    Ac(6,7) =  g * versinc(phi);
    Ac(6,8) =  g * versinc(theta) * cos(phi);

    % Euler attitude kinematics: euler_dot = J(phi,theta) * omega
    J = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
    Ac(7:9,10:12) = J;

    % Gyroscopic coupling (Euler's rigid-body equations, SDC form)
    Ac(10,11) = -x(12) * (qp.I_zz - qp.I_yy) / qp.I_xx;
    Ac(11,12) = -x(10) * (qp.I_xx - qp.I_zz) / qp.I_yy;
    Ac(12,10) = -x(11) * (qp.I_yy - qp.I_xx) / qp.I_zz;
end

function x = versinc(x)
%               1 - cos(x)
% versinc(x) := ----------
%                   x
    x = sin(x/2) * sinc1(x/2);
end