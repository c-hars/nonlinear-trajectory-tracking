

function [r_,x_ref_fcn,tgrid,tspan,x0] = load_fig8_traj(qp,C,opts)

    arguments
        qp
        C
        opts.t_maneuver = 9.0
        opts.AttitudeRepresentation = 'Euler';
    end

    t_maneuver = ceil(opts.t_maneuver/qp.Ts)*qp.Ts;
    w = 2*pi/t_maneuver;
    x_ref_fcn = @(t) 4*[sin(w*t); -cos(w*t)*sin(w*t); cos(w*t); ...
                      w*cos(w*t); -w*cos(2*w*t); -w*sin(w*t); ...
                      zeros(6,1)];
    % w = 2.0 * 2*pi/t_maneuver;
    % x_ref_fcn = @(t) 4*[sin(w*t); -1.5*cos(w*t)*sin(w*t); -cos(0.5*w*t); ...
    %                   w*cos(w*t); -1.5*w*cos(2*w*t); w*sin(0.5*w*t); ...
    %                   zeros(6,1)];

    x0 = x_ref_fcn(0);
    % x0(4:6) = zeros(3,1); % enable this -> copter starts at correct position, but velocity mismatch with trajectory

    switch lower(opts.AttitudeRepresentation)
        case 'quaternion'
            x0 = [x0(1:6); [1;0;0;0]; zeros(3,1)];
        otherwise
            x0 = [x0(1:6); zeros(6,1)];
    end

    tgrid = 0:qp.Ts:t_maneuver;
    r_ = arrayfun(x_ref_fcn, tgrid, 'UniformOutput', false);
    r_ = cell2mat(r_);
    r_ = C * r_;

    tspan = [0 t_maneuver];

end