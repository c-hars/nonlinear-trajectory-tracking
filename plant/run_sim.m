function [T,X,U,compute_time_ctrl,compute_time_ode45,solve_info,U_raw,status] = run_sim(qp, tspan, x0, enabled_fcn, ctrl_fcn, opts)
    arguments
        qp
        tspan
        x0
        enabled_fcn
        ctrl_fcn = @(t,x,k) qp.K * x
        opts.AttitudeRepresentation = 'euler'
        opts.MaxWallClock = inf     % [s] per-run compute budget
        opts.Integrator = []        % 'ode45' or 'rk4'
    end

    if isempty(opts.Integrator)
        if qp.Ts < 1/20
            opts.Integrator = 'rk4';
        else
            opts.Integrator = 'ode45';
        end
    else
        if qp.Ts < 1/20 && strcmpi(opts.Integrator,'rk4')
            warning('Overriding rk4 request with ode45 (rk4 is inaccurate below 20Hz loop rate: current loop rate is %gHz)', 1/qp.Ts)
            opts.Integrator = 'ode45';
        end
    end

    M  = 1;  % number of sub-samples per control interval (M >= 1)
    dt = qp.Ts;
    N  = ceil(tspan(end) / dt);
    nx = numel(x0);
    nu = qp.n_rotors;

    T = zeros(N*M + 1, 1);
    X = zeros(N*M + 1, nx);
    U = zeros(N, nu);
    U_raw = zeros(N, nu);
    compute_time_ctrl  = zeros(N, 1);
    compute_time_ode45 = zeros(N, 1);
    solve_info = cell(N, 1);
    T(1) = 0;
    X(1,:) = x0;

    status = struct('diverged', false, 'k_abort', NaN, 't_abort', NaN, 'message', "");

    switch lower(opts.AttitudeRepresentation)
        case 'quaternion'
            assert(nx == 13);
        otherwise
            assert(nx == 12);
    end

    x = x0;
    run_clock = tic;
    for k = 1:N

        if toc(run_clock) > opts.MaxWallClock
            status.diverged    = true;
            status.k_abort     = k;
            status.t_abort     = t;
            status.message     = sprintf("Wall-clock budget %.0fs exceeded at k=%d", opts.MaxWallClock, k);
            status.identifier  = "ctrl:wallclock";
            warning("ABORTING: %s", status.message)
            n_keep = (k-2)*M + 1;
            [T,X,U,U_raw,compute_time_ctrl,compute_time_ode45,solve_info] = ...
                truncate_data(n_keep, k, T, X, U, U_raw, compute_time_ctrl, compute_time_ode45, solve_info);
            return
        end

        t = (k-1) * dt;

        % Compute control input (ZOH for next Ts)
        clock_start = tic;
        if k == 1, u_prev = zeros(nu,1); else, u_prev = U(k-1,:)'; end

        args = {t, x, k};
        if nargin(ctrl_fcn) == 4, args{end+1} = u_prev; end

        try
            [u, solve_info{k}] = ctrl_fcn(args{:});
        catch err
            warning(err.identifier, "Error in ctrl_fcn: %s", err.message)
            status.diverged    = true;
            status.k_abort     = k;
            status.t_abort     = t;
            status.message     = string(err.message);
            status.identifier  = string(err.identifier);
            n_keep = (k-1)*M + 1;
            [T,X,U,U_raw,compute_time_ctrl,compute_time_ode45,solve_info] = ...
                truncate_data(n_keep, k, T, X, U, U_raw, compute_time_ctrl, compute_time_ode45, solve_info);
            return
        end
        compute_time_ctrl(k) = toc(clock_start);

        % Apply actuator limits
        U_raw(k,:) = u;
        u = min(u, qp.max_du(:));
        u = max(u, -qp.nominal_omegas(:));  % floor at 0 RPM (positive thrust only)
        U(k,:) = u;

        
        % Integrate from t to t+Ts, u held constant (ZOH)
        clock_start = tic;
        qp.enabled = enabled_fcn(t);
        dynamics = @(t,x) nonlinear_dynamics(x, u, qp, opts.AttitudeRepresentation);
        switch opts.Integrator
            case 'rk4'
                x_sub = x;
                dt_sub = dt / M;
                for m = 1:M
                    t_sub = t + (m-1)*dt_sub;
                    x_sub = rk4_step(dynamics, t_sub, x_sub, dt_sub);
                    idx = (k-1)*M + 1 + m;
                    T(idx)   = t_sub + dt_sub;
                    X(idx,:) = x_sub';
                end
                x = x_sub;

            case 'ode45'
                t_query = linspace(t, t+dt, M+1);
                [~, x_local] = ode45(dynamics, t_query, x);
                if M == 1
                    x_local = x_local([1 end],:);
                end
                idx = (k-1)*M + 2 : k*M + 1;
                T(idx)   = t_query(2:end)';
                X(idx,:) = x_local(2:end,:);
                x = x_local(end,:)';

        end
        compute_time_ode45(k) = toc(clock_start);

    end

    solve_info = [solve_info{:}];

end


function [T_out,X_out,U_out,U_raw_out,compute_time_ctrl_out,compute_time_ode45_out,solve_info_out] = ...
        truncate_data(n_keep, k, T, X, U, U_raw, compute_time_ctrl, compute_time_ode45, solve_info)
    T_out                  = T(1:n_keep);
    X_out                  = X(1:n_keep, :);
    U_out                  = U(1:k-1, :);
    U_raw_out              = U_raw(1:k-1, :);
    compute_time_ctrl_out  = compute_time_ctrl(1:k-1);
    compute_time_ode45_out = compute_time_ode45(1:k-1);
    solve_info_out         = solve_info(1:k-1);
    solve_info_out         = [solve_info_out{:}];
end

function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end