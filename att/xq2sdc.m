function xmap = xq2sdc(rep)
% Takes xq, 13-element quaternion plant state -> 12-state SDC layout  [pos; vel; att(3); omega] (SDC-controller-native coordinates - att(3) is either Euler/Quat_vec/MRP/FRA).
% Example usage: ctrl_fcn = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, ...)
    switch lower(rep)
        case 'euler'
            xmap = @(x) [x(1:6); quat2euler_(x(7:10)); x(11:13)];
        case 'mrp'
            xmap = @(x) [x(1:6); quat2mrp_(x(7:10));   x(11:13)];
        case 'fra'
            xmap = @(x) [x(1:6); quat2fra_(x(7:10));   x(11:13)];
        case 'quaternion'
            xmap = @(x) [x(1:6); quat2qvec_(x(7:10));  x(11:13)];
        otherwise
            error("xq2sdc: unknown representation '%s'", rep);
    end
end

function q = qfix(q)
% Normalise + principal branch (q0 >= 0). Shared by all converters.
    q = q / norm(q);
    if q(1) < 0, q = -q; end
end

function qv = quat2qvec_(q)
% Vector part only; SDC recovers q0 = +sqrt(1 - |qv|^2) internally, hence the q0 >= 0 branch fix here is required (not cosmetic).
    q  = qfix(q);
    qv = q(2:4);
end

function p = quat2mrp_(q)
% Principal MRP set (|p| <= 1 guaranteed by q0 >= 0).
    q = qfix(q);
    p = q(2:4) / (1 + q(1));
end

function th = quat2fra_(q)
% Principal rotation vector (angle <= pi).
% NB: NOT branch-continuous.
    q  = qfix(q);
    nv = norm(q(2:4));
    if nv < 1e-12
        th = zeros(3,1);
    else
        th = (2*atan2(nv, q(1)) / nv) * q(2:4);
    end
end

function e = quat2euler_(q)
    q  = qfix(q);
    b0 = q(1); bv = q(2:4);
    C_bI = (b0^2 - bv.'*bv)*eye(3) + 2*(bv*bv.') - 2*b0*S(bv);
    phi   = atan2(C_bI(2,3), C_bI(3,3));
    theta = atan2(-C_bI(1,3), sqrt(C_bI(1,1)^2 + C_bI(1,2)^2));
    psi   = atan2(C_bI(1,2), C_bI(1,1));
    e = [phi; theta; psi];
end