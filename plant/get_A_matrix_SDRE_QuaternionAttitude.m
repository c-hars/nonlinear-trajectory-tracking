function A = get_A_matrix_SDRE_QuaternionAttitude(x,qp)

    % State: [pos(1:3); vel(4:6); q1,q2,q3(7:9); wx,wy,wz(10:12)]
    % q0 is recovered implicitly
    
    q1 = x(7);  q2 = x(8);  q3 = x(9);
    q0 = sqrt(1 - q1^2 - q2^2 - q3^2); % Recover q0 from norm constraint
    persistent haswarned
    if isempty(haswarned) && q0 < 0.1
        warning("q0 dipped below 1e-1");
        haswarned=true;
    end

    wx = x(10); wy = x(11); wz = x(12);

    g  = 9.81;

    A = zeros(12);

    % --- Position kinematics ---
    C_Ib = [(q0^2+q1^2-q2^2-q3^2),      2*(q1*q2 - q0*q3),      2*(q1*q3 + q0*q2);
                2*(q1*q2 + q0*q3),  (q0^2-q1^2+q2^2-q3^2),      2*(q2*q3 - q0*q1);
                2*(q1*q3 - q0*q2),      2*(q2*q3 + q0*q1),  (q0^2-q1^2-q2^2+q3^2)];
    A(1:3,4:6) = C_Ib;
    

    % --- Coriolis ---
    % Best left out (same treatment as the Euler formulation)
    % A(4:6, 4:6) = -S([wx;wy;wz]);

    % --- Gravity residual ---
    % Hover: simplifies/linearises to:
    %              ( 2*g)*q2
    %              (-2*g)*q1
    %                      0
    % so SDC form should resemble this near q=[1 0 0 0]
    % Nonlinear:
    %               g*(2*q0*q2 - 2*q1*q3)
    %              -g*(2*q0*q1 + 2*q2*q3)
    %               g*(1 - q0^2 + q1^2 + q2^2 - q3^2)
    % Neglect the affine term (1-q0^2)? Approximately 0, and is an affine term - can't be pulled into the SDC form too readily.
    % Then we get the SDC factorisation:
    % A(4:6, 7:9) =  g * [-2*q3,  2*q0,   0;
    %                     -2*q0, -2*q3,   0;
    %                        q1,    q2, -q3];    
    % 
    % ... but q0^2 = 1 - q1^2 - q2^2 - q3^2. So 1-q0^2 = q1^2 + q2^2 + q3^2,
    % which is SDC factorisable, and the whole thing becomes quite simple:
    % A(4:6, 7:9) =  g * [-2*q3,  2*q0,   0;
    %                     -2*q0, -2*q3,   0;
    %                      2*q1,  2*q2,   0];
    % Also try this "symmetric split" - slightly better behaved in general
    A(4:6, 7:9) =  g * [  -q3, 2*q0, -q1;
                        -2*q0,  -q3, -q2;
                         2*q1, 2*q2,   0];
    % A(4:6, 7:9) =  g * [ 0, 2, 0;
    %                     -2, 0, 0;
    %                      0, 0, 0 ];
    % this linearish approximation actually does work well tho, better in some respects, makes xy tracking better, z tracking worse?
    

    % --- Quaternion kinematics ---
    % A(7:9, 10:12) = 0.5 * [ q0, -q3,  q2;    % dq1/dt
    %                         q3,  q0, -q1;    % dq2/dt
    %                        -q2,  q1,  q0];   % dq3/dt

    % NB: below yields slightly tracking that is (slightly) better but more aggressive
    % The simple version above is better for low sample rates (20Hz), capable of handling full range of aggressive maneuvers (up to 3.6 sin amplitude despite 20Hz), just with slightly less performant tracking.
    % The below version is better for eeking performance out - pretty noticeable tbf - and worth it if the maneuvers aren't aggressive, or the sample rate is high enough (then aggressive maneuvers can be handled just fine, very well actually, better than above)
    Omega_w = [   0,  wz, -wy;
                -wz,   0,  wx;
                 wy, -wx,   0 ];
    Xi_q    = [  q0, -q3,  q2;
                 q3,  q0, -q1;
                -q2,  q1,  q0 ];
    % alpha = 0.99; % Convex blend (alpha in [0,1])
    alpha = 1; % Convex blend (alpha in [0,1])
    A(7:9, 7:9)   = alpha     * 0.5*Omega_w;
    A(7:9, 10:12) = (1-alpha) * 0.5*Xi_q + alpha * 0.5*q0*eye(3);
    % A(7:9, 10:12) = eye(3);

    % --- Gyroscopic term (angular velocity dynamics) ---
    % Usually trivial and can be left out
    A(10,11) = -x(12)*(qp.I_zz-qp.I_yy)/qp.I_xx;
    A(11,12) = -x(10)*(qp.I_xx-qp.I_zz)/qp.I_yy;
    A(12,10) = -x(11)*(qp.I_yy-qp.I_xx)/qp.I_zz;

end