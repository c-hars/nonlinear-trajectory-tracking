function dxdt = nonlinear_dynamics(x, u, qp, opts)
    arguments
        x (:,1)
        u (:,1)
        qp struct
        opts.AttitudeRepresentation = 'euler'
    end

    if ~isreal(u), u = real(u); warning("u was not real"); end
    if ~isreal(x), x = real(x); warning("x was not real"); end

    % Motor mixing (shared across all branches)
    M_mix_NL = [ qp.kF * qp.enabled;
                 qp.kF * qp.enabled .* qp.y;
                -qp.kF * qp.enabled .* qp.x;
                -qp.kM * qp.enabled .* qp.dirs ];

    abs_angvels = u + qp.nominal_omegas';
    abs_angvels = max(abs_angvels, 0);
    abs_angvels = min(abs_angvels, qp.max_RPM * 2*pi/60);

    F_t = M_mix_NL(1,:) * abs_angvels.^2;
    T_c = M_mix_NL(2:4,:) * abs_angvels.^2;

    g = 9.81;
    I = diag([qp.I_xx, qp.I_yy, qp.I_zz]);

    if strcmpi(opts.AttitudeRepresentation, 'euler')

        v_c   = x(4:6);
        phi   = x(7); theta = x(8); psi = x(9);
        omega = x(10:12);

        C_bI = C_x(phi) * C_y(theta) * C_z(psi);
        C_Ib = C_bI';

        v_c_dot = -S(omega)*v_c + (1/qp.m)*([0;0;F_t] + C_bI*[0;0;-qp.m*g]);

        J = [1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
             0, cos(phi),            -sin(phi);
             0, sin(phi)*sec(theta),  cos(phi)*sec(theta)];
        euler_dot = J * omega;

        omega_dot = -I\(S(omega)*I*omega) + I\T_c;

        dxdt = zeros(12,1);
        dxdt(1:3)   = C_Ib * v_c;
        dxdt(4:6)   = v_c_dot;
        dxdt(7:9)   = euler_dot;
        dxdt(10:12) = omega_dot;

    elseif strcmpi(opts.AttitudeRepresentation, 'quaternion')
        % State: [pos(1:3); vel_b(4:6); q0,q1,q2,q3(7:10); omega(11:13)]

        v_c   = x(4:6);
        q0    = x(7); q1 = x(8); q2 = x(9); q3 = x(10);
        omega = x(11:13);

        % Renormalise defensively
        nrm = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
        q0 = q0/nrm; q1 = q1/nrm; q2 = q2/nrm; q3 = q3/nrm;

        C_Ib = [(q0^2+q1^2-q2^2-q3^2),  2*(q1*q2 - q0*q3),      2*(q1*q3 + q0*q2);
                 2*(q1*q2 + q0*q3),     (q0^2-q1^2+q2^2-q3^2),  2*(q2*q3 - q0*q1);
                 2*(q1*q3 - q0*q2),      2*(q2*q3 + q0*q1),     (q0^2-q1^2-q2^2+q3^2)];
        C_bI = C_Ib';

        v_c_dot = -S(omega)*v_c + (1/qp.m)*([0;0;F_t] + C_bI*[0;0;-qp.m*g]);

        % Quaternion kinematics: q_dot = (1/2) * Xi(q) * [0; omega]
        Xi = [q0, -q1, -q2, -q3;
              q1,  q0, -q3,  q2;
              q2,  q3,  q0, -q1;
              q3, -q2,  q1,  q0];
        q_dot = 0.5 * Xi * [0; omega];

        omega_dot = -I\(S(omega)*I*omega) + I\T_c;

        dxdt = zeros(13,1);
        dxdt(1:3)   = C_Ib * v_c;
        dxdt(4:6)   = v_c_dot;
        dxdt(7:10)  = q_dot;
        dxdt(11:13) = omega_dot;

    elseif strcmpi(opts.AttitudeRepresentation, 'mrp')
        % State: [pos(1:3); vel_b(4:6); p1,p2,p3(7:9); omega(10:12)]

        v_c   = x(4:6);
        p     = x(7:9);
        omega = x(10:12);

        p_sq = p'*p;
        Sp   = S(p);
        C_Ib = eye(3) + (4*(1 - p_sq)/(1+p_sq)^2)*Sp + (8/(1+p_sq)^2)*(Sp*Sp);
        C_bI = C_Ib';

        v_c_dot = -S(omega)*v_c + (1/qp.m)*([0;0;F_t] + C_bI*[0;0;-qp.m*g]);

        % MRP kinematics: p_dot = G(p) * omega
        G = 0.5 * ((1 - p_sq)/2 * eye(3) + Sp + p*p');
        p_dot = G * omega;

        omega_dot = -I\(S(omega)*I*omega) + I\T_c;

        dxdt = zeros(12,1);
        dxdt(1:3)   = C_Ib * v_c;
        dxdt(4:6)   = v_c_dot;
        dxdt(7:9)   = p_dot;
        dxdt(10:12) = omega_dot;

    elseif strcmpi(opts.AttitudeRepresentation, 'fra')
        % State: [pos(1:3); vel_b(4:6); th1,th2,th3(7:9); omega(10:12)]

        v_c   = x(4:6);
        th    = x(7:9);
        omega = x(10:12);

        Sth = S(th);
        ang = norm(th);
        g_  = sinc1(ang);
        f_  = 0.5 * sinc1(ang/2)^2;
        if ang < 1e-2
            h_ = (1/6) - ang^2/120 + ang^4/5040;
            % h_ = ((1/6) - (11/2250)*ang^2) / (1 + (1/42)*ang^2);
        else
            h_ = (ang - sin(ang)) / ang^3;
        end

        C_Ib = eye(3) + g_*Sth + f_*(Sth*Sth);
        C_bI = C_Ib';

        v_c_dot = -S(omega)*v_c + (1/qp.m)*([0;0;F_t] + C_bI*[0;0;-qp.m*g]);

        % FRA kinematics: theta_dot = inv(F) * omega
        F = eye(3) - f_*Sth + h_*(Sth*Sth);
        fra_dot = F \ omega;

        omega_dot = -I\(S(omega)*I*omega) + I\T_c;

        dxdt = zeros(12,1);
        dxdt(1:3)   = C_Ib * v_c;
        dxdt(4:6)   = v_c_dot;
        dxdt(7:9)   = fra_dot;
        dxdt(10:12) = omega_dot;
    end
end