
clear
load_paths
load_copter_params

%% Main parameters.

% Main parameters to experiment with are here.
% You can also adjust the cost matrices - please do so in the get_weights() function ("plant/get_weights.m").

qp.Ts = 1/20;         % sample rate. NB: qp stands for "quadcopter parameters" (the plant was originally a quadcopter:))
maneuver_time = 10.0;  % time the maneuever needs to be completed in [seconds]

SDCAttRep = 'Euler';    % for SDOPT. Choose from: Euler, Quaternion, MRP, FRA.
SDOPTPreviewHorizon = 2.0; % [seconds]

%% Linear control (Linear Quadratic Tracking)

Ac      = get_A_matrix();
Bc      = get_B_matrix(qp);
[Ad,Bd] = c2d_zoh_expm(Ac, Bc, qp.Ts);
[Q, R, C, Qy, Qyf, ~, sel] = get_weights(qp, 'Euler');
[r_,x_ref_fcn,~,tspan,x0]  = load_fig8_traj(qp, C, ...
    t_maneuver=maneuver_time, AttitudeRepresentation='Euler');

t_start = tic;
[K_, xref_, ff_] = solve_LQT(r_, length(r_), Ad, Bd, C, Qy, Qyf, R);
lqt_precompute = toc(t_start);

ctrl_fcn   = @(t,x,k) deal( ...
    K_(:,:,k) * (x - xref_(:,k)) + ff_(:,k), ...
    struct('ExitFlag', 1));
thrust_fcn = @(t) ones(1,6);

[t, X, U, c1] = run_sim(qp, tspan, x0, thrust_fcn, ctrl_fcn, ...
    AttitudeRepresentation='Euler');

[J, Jy, Ju, Jy_i] = compute_J_LQT_v3(t, X, U, Qy, Qyf, R, C, x_ref_fcn, qp.Ts, ...
    AttitudeRepresentation='Euler');

fprintf('\n--- LQT ---\n')
fprintf('  Compute : %.2f s (precompute) + %.2f s (sim overhead)\n', lqt_precompute, sum(c1))
print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)

figure(1); clf
do_plots(t, X, U, r_, qp)
sgtitle('LQT')


%% Nonlinear control (SD-OPT, "State Dependent Optimal Preview Tracking")

% Regardless of the controller, the plant always uses quaternion dynamics under the hood
% xq2sdc maps the 13-state quaternion plant to controller-native coordinates
% get_weights() recomputes Qy, Qyf to ensure cost-function equivalence across attitude representations (and sample rates)
xmap = xq2sdc(SDCAttRep);
[Q, R, C, Qy, Qyf, ~, sel] = get_weights(qp, SDCAttRep);
[r_, x_ref_fcn, tgrid, tspan, x0] = load_fig8_traj(qp, C, ...
    t_maneuver=maneuver_time, AttitudeRepresentation=SDCAttRep);

A_fcn = get_SDC_A_function(SDCAttRep);

% SD-OPT control function
ctrl_fcn = @(t,x,k,u_prev) compute_u_SDOPT(t, xmap(x), k, u_prev, ...
    r_, C, Qy, R, Qyf, qp, ...
    PreviewHorizon = SDOPTPreviewHorizon, ...
    SDC_A_function = A_fcn, ...
    SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

thrust_fcn = @(t) ones(1,6);
x0_q = [x0(1:6); [1;0;0;0]; zeros(3,1)];

[t, X, U, c1, ~, ~, U_raw] = run_sim(qp, tspan, x0_q, thrust_fcn, ctrl_fcn, ...
    AttitudeRepresentation='Quaternion');

% Evaluate cost in Euler coordinates regardless of controller representation
[~, ~, C_eul, Qy_eul, Qyf_eul] = get_weights(qp, 'Euler');
[J, Jy, Ju, Jy_i] = compute_J_LQT_v3(t, X, U, Qy_eul, Qyf_eul, R, C_eul, ...
    x_ref_fcn, qp.Ts, AttitudeRepresentation='Quaternion');

fprintf('\n--- SD-OPT (with %s attitude) ---\n', SDCAttRep)
fprintf('  Compute : %.2f s total → %.1f Hz equivalent\n', sum(c1), length(U)/sum(c1))
print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)
print_saturation(U_raw, qp)

figure(2); clf
do_plots(t, X, U, r_, qp)
sgtitle(sprintf('SD-OPT (with %s attitude)', SDCAttRep))


%% Compare with precomputed NLMPC (too slow to recompute live)

nlmpc_file = sprintf('data/nlmpc_fra_cold_con_Ts%gHz_tman%s_H%g.mat', 1/qp.Ts, ...
    strrep(num2str(maneuver_time), '.', 'p'), SDOPTPreviewHorizon);

% if qp.Ts ~= 1/20 || SDOPTPreviewHorizon ~= 2.0
%     fprintf('\n--- NLMPC: precomputed results are at 20 Hz with Horizon 2.0s, current Ts is 1/%g Hz and Horizon %gs. Skipping ---\n', 1/qp.Ts, SDOPTPreviewHorizon)
if isfile(nlmpc_file)
    nlmpc_data = load(nlmpc_file);

    [~, ~, C_eul, Qy_eul, Qyf_eul] = get_weights(qp, 'Euler');
    [J, Jy, Ju, Jy_i] = compute_J_LQT_v3(nlmpc_data.t, nlmpc_data.X, nlmpc_data.U, ...
        Qy_eul, Qyf_eul, R, C_eul, x_ref_fcn, qp.Ts, ...
        AttitudeRepresentation='Quaternion');
    fprintf('\n--- NLMPC (default Q/R matrices, FRA attitude, cold start, actuator constraints) ---\n')
    fprintf('  Compute : %.2f s (total), %.2f s (first iter) + %.2fs/%.2fs/%.2fs/%.2fs (iter-wise mean/median/min/max thereafter)\n', ...
        sum(nlmpc_data.compute_time_ctrl(1:end)), ...
        nlmpc_data.compute_time_ctrl(1), ...
        mean(nlmpc_data.compute_time_ctrl(2:end)), ...
        median(nlmpc_data.compute_time_ctrl(2:end)), ...
        min(nlmpc_data.compute_time_ctrl(2:end)), ...
        max(nlmpc_data.compute_time_ctrl(2:end)))
    print_tracking_metrics(nlmpc_data.X, r_, J, Jy, Ju, Jy_i)
    print_saturation(nlmpc_data.U, qp, 0.005)

    figure(3); clf
    do_plots(nlmpc_data.t, nlmpc_data.X, nlmpc_data.U, r_, qp)
    sgtitle(sprintf('NLMPC (with FRA attitude)'))
else
    fprintf('\n--- NLMPC: no precomputed result for Ts=%gHz, Horizon=%gs, t_man=%.2f s ---\n', 1/qp.Ts, SDOPTPreviewHorizon, maneuver_time)
    figure(3); clf
end


%% ---- Local functions ----

function print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)
    N = length(r_);
    e_pos  = X(1:N, 1:3) - r_(1:3, 1:N)';
    rmse   = 100 * sqrt(mean(e_pos.^2, 'all'));
    e_term = 100 * norm(e_pos(end,:));

    Jy_pos = sum(Jy_i(1:3));

    fprintf('  Cost    : %.2f  (tracking: %.2f, input: %.2f)\n', J, Jy, Ju)
    fprintf('    of which position (xyz): %.2f  (%.0f%% of tracking cost, %.0f%% overall)\n', Jy_pos, 100*Jy_pos/Jy, 100*Jy_pos/J)
    % fprintf('    per-state contributions to tracking cost: [ %.0f, %.0f, %.0f, %.0f, %.0f, %.0f ] %% for states [ x, y, z, omega_3, omega_1_dot, omega_2_dot ]\n', 100*Jy_i/sum(Jy_i))
    fprintf('  RMSE    : %.2f cm\n', rmse)
    fprintf('  Terminal: %.2f cm\n', e_term)
end

function print_saturation(U_raw, qp, tol)
    if nargin < 3, tol = 0; end
    upper =  qp.max_du    * (1 - tol);
    lower = -qp.nominal_omegas * (1 - tol);
    sat = (U_raw > upper) | (U_raw < lower);
    if tol > 0
        fprintf('  Saturation (within %.1f%% of bounds): %.1f%% of actuator samples, %.1f%% of timesteps\n', ...
            100*tol, 100*sum(sat,'all')/numel(sat), 100*sum(any(sat,2))/size(sat,1))
    else
        fprintf('  Saturation: %.1f%% of actuator samples, %.1f%% of timesteps\n', ...
            100*sum(sat,'all')/numel(sat), 100*sum(any(sat,2))/size(sat,1))
    end
end

function A_fcn = get_SDC_A_function(rep)
    switch lower(rep)
        case 'euler',      A_fcn = @(xk,qp) get_A_matrix_SDRE_EulerAttitude(xk,qp);
        case 'mrp',        A_fcn = @(xk,qp) get_A_matrix_SDRE_MRPAttitude(xk,qp);
        case 'fra',        A_fcn = @(xk,qp) get_A_matrix_SDRE_FRAAttitude(xk,qp);
        case 'quaternion', A_fcn = @(xk,qp) get_A_matrix_SDRE_QuaternionAttitude(xk,qp);
    end
end