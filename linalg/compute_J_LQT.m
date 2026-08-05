function [J,Jy,Ju,Ji] = compute_J_LQT(t, x_, u_, Qy, Qyf, R, C, x_ref_fcn, Ts, opts)
% Evaluate LQT cost in a common metric across attitude representations.
% Each state sample is converted to the Euler-equivalent 12-state layout
% [p; v_c; euler; omega] before C is applied, so J is directly comparable
% across runs regardless of representation.
%
% Output errors on pure-angle channels (C-row support only in cols 7:9)
% are wrapped to (-pi, pi]. Set WrapAngleErrors=false to disable.

    arguments
        t; x_; u_; Qy; Qyf; R; C; x_ref_fcn; Ts
        opts.AttitudeRepresentation (1,1) string ...
            {mustBeMember(opts.AttitudeRepresentation, ...
             ["Euler","FRA","MRP","Quaternion"])} = "Euler"
        opts.WrapAngleErrors (1,1) logical = true
    end
    rep = opts.AttitudeRepresentation;

    N  = size(u_,1);
    ny = size(C,1);

    tk_ = ((0:N)') * Ts;
    xk_ = interp1(t, x_, min(tk_, t(end)));

    if ~(norm(Qy - diag(diag(Qy))) < 1e-6)
        warning(['Ji breakdown assumes diagonal Qy; J, Jy, Ju are ' ...
                 'still accurate.'])
    end

    % Output channels eligible for angle wrapping: C-row support only in
    % the Euler-angle columns of the Euler-equivalent state.
    ang_cols  = 7:9;
    oth_cols  = setdiff(1:12, ang_cols);
    wrap_rows = opts.WrapAngleErrors & ...
                any(C(:,ang_cols) ~= 0, 2) & ~any(C(:,oth_cols) ~= 0, 2);

    Jy = 0;
    Ju = 0;
    Ji = zeros(ny,1);
    Qy_diag = diag(Qy);

    for k = 1:N
        rk = C * x_ref_fcn(tk_(k));
        xe = to_euler_state(xk_(k,:)', rep);
        ek = wrap_err(C*xe - rk, wrap_rows);
        uk = u_(k,:)';

        Jy = Jy + ek' * Qy * ek;
        Ju = Ju + uk' * R * uk;
        Ji = Ji + ek.^2 .* Qy_diag;
    end

    % Terminal cost
    xe = to_euler_state(xk_(N+1,:)', rep);
    ef = wrap_err(C*xe - C*x_ref_fcn(tk_(N+1)), wrap_rows);
    Jy = Jy + ef' * Qyf * ef;
    Ji = Ji + ef.^2 .* diag(Qyf);

    J = Jy + Ju;
end


function e = wrap_err(e, wrap_rows)
    e(wrap_rows) = mod(e(wrap_rows) + pi, 2*pi) - pi;
end