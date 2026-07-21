function e = quat2euler(q)
% Scalar-first quaternion -> 3-2-1 Euler via C_bI entries.
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    C11 = q0^2 + q1^2 - q2^2 - q3^2;
    C12 = 2*(q1*q2 + q0*q3);
    C13 = 2*(q1*q3 - q0*q2);
    C23 = 2*(q2*q3 + q0*q1);
    C33 = q0^2 - q1^2 - q2^2 + q3^2;
    e = [ atan2(C23, C33);
          atan2(-C13, sqrt(C11^2 + C12^2));
          atan2(C12, C11) ];
end