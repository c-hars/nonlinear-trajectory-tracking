
clear
load_paths
load_copter_params

%% Main parameters.

% Main parameters to experiment with are here.
% Other parameters can also be changed - e.g. modify the trajectory at ("utils/load_fig_8.m") and cost matrices at ("plant/get_weights.m").

qp.Ts = 1/20;  % sample rate. NB: qp stands for "quadcopter parameters" (the plant was originally a quadcopter).
maneuver_time = 10.0;  % time the maneuever needs to be completed in [seconds]

SDCAttRep = 'FRA'; % for SDOPT. Choose from: Euler, Quaternion, MRP, FRA.
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

[J, Jy, Ju, Jy_i] = compute_J_LQT(t, X, U, Qy, Qyf, R, C, x_ref_fcn, qp.Ts, ...
    AttitudeRepresentation='Euler');

fprintf('\n--- LQT ---\n')
fprintf('  Compute : %.2f s (precompute) + %.2f s (sim overhead)\n', lqt_precompute, sum(c1))
print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)

figure(1); clf
do_plots(t, X, U, r_, qp, 'Euler')
sgtitle('LQT')


%% Nonlinear control (SD-OPT, "State Dependent Optimal Preview Tracking")

% Regardless of the controller, the plant always uses quaternion dynamics under the hood.
% xq2sdc maps the 13-state quaternion plant to controller-native coordinates.
% get_weights() recomputes Qy, Qyf to ensure cost-function equivalence across attitude representations (and sample rates).


SimDataRep = 'Quaternion';
xmap = xq2sdc(SDCAttRep);

[Q, R, C, Qy, Qyf, ~, sel] = get_weights(qp, SDCAttRep);
[r_, x_ref_fcn, tgrid, tspan, x0] = load_fig8_traj(qp, C, ...
    t_maneuver=maneuver_time, AttitudeRepresentation=SimDataRep);

A_fcn = get_SDC_A_function(SDCAttRep);

% SD-OPT control function
ctrl_fcn = @(t,x,k,u_prev) compute_u_SDOPT(t, xmap(x), k, u_prev, ...
    r_, C, Qy, R, Qyf, qp, ...
    PreviewHorizon = SDOPTPreviewHorizon, ...
    SDC_A_function = A_fcn, ...
    SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

thrust_fcn = @(t) ones(1,6);
[t, X, U, c1, c2, stats, U_raw] = run_sim(qp, tspan, x0, thrust_fcn, ctrl_fcn, ...
    AttitudeRepresentation=SimDataRep);

% Evaluate cost in Euler coordinates regardless of controller representation
[~, ~, C_eul, Qy_eul, Qyf_eul] = get_weights(qp, 'Euler');
[J, Jy, Ju, Jy_i] = compute_J_LQT(t, X, U, Qy_eul, Qyf_eul, R, C_eul, ...
    x_ref_fcn, qp.Ts, AttitudeRepresentation=SimDataRep);

fprintf('\n--- SD-OPT (with %s attitude) ---\n', SDCAttRep)
fprintf('  Compute : %.2f s total → %.1f Hz equivalent\n', sum(c1), length(U)/sum(c1))
fprintf('    running the simulation (rk4 or ode45) took up the other %.2f s\n', sum(c2))
print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)
print_saturation(U_raw, qp)

figure(2); clf
do_plots(t, X, U, r_, qp, SimDataRep)
sgtitle(sprintf('SD-OPT (with %s attitude)', SDCAttRep))


%% Alternative: 3D trajectory plot - show the reference path versus and what was actually tracked.

% n_disp = round(maneuver_time * 50);

% nr = size(r_,2);
% t = t(1:nr);
% X = X(1:nr,:);

% ri = interp1(t,r_(1:3,:)',linspace(0,t(end),n_disp));
% xi = interp1(t,X(:,1:3),linspace(0,t(end),n_disp));

% figure(2); clf
% plot3(ri(:,1),ri(:,2),ri(:,3)); axis equal; grid on
% hold on; scatter3(xi(:,1), xi(:,2), xi(:,3), 4, parula(size(xi,1)), 'o')
% xlabel('X'); ylabel('Y'); zlabel('Z')
% view(125,40)


%% Compare with precomputed NLMPC (too slow to recompute live)

nlmpc_file = sprintf('data/nlmpc_fra_cold_con_Ts%gHz_tman%s_H%g.mat', 1/qp.Ts, ...
    strrep(num2str(maneuver_time), '.', 'p'), SDOPTPreviewHorizon);

if isfile(nlmpc_file)
    nlmpc_data = load(nlmpc_file);
    SimDataRep = 'Quaternion';

    [~, ~, C_eul, Qy_eul, Qyf_eul] = get_weights(qp, 'Euler');
    [J, Jy, Ju, Jy_i] = compute_J_LQT(nlmpc_data.t, nlmpc_data.X, nlmpc_data.U, ...
        Qy_eul, Qyf_eul, R, C_eul, x_ref_fcn, qp.Ts, ...
        AttitudeRepresentation=SimDataRep);
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
    do_plots(nlmpc_data.t, nlmpc_data.X, nlmpc_data.U, r_, qp, SimDataRep)
    sgtitle(sprintf('NLMPC (with FRA attitude)'))
else
    fprintf('\n--- NLMPC: no precomputed result for Ts=%gHz, Horizon=%gs, t_man=%.2f s ---\n', 1/qp.Ts, SDOPTPreviewHorizon, maneuver_time)
    figure(3); clf % clear any old/invalid plot
end


%% SDOPT w/ high loop rates

qp.Ts = 1/1000;
maneuver_time = 5.0;
SDCAttRep = 'Quaternion';

% Decimation factors
% df = 2.^(0:1:floor(log2((1/20)/qp.Ts)));  % geometric progression, up to ~20Hz preview
df = 2.^[0, floor(log2((1/20)/qp.Ts))];  % straight progression: fine-rate then coarse-rate (no intermediate geometric progression)
dopts = struct('DecimationFactors', df, 'StepsPerTier', df(end));

% Problem setup
xmap = xq2sdc(SDCAttRep);
A_fcn = get_SDC_A_function(SDCAttRep);
[Q, R, C, Qy, Qyf, ~, sel] = get_weights(qp, SDCAttRep);
[r_, x_ref_fcn, tgrid, tspan, x0] = load_fig8_traj(qp, C, ...
    t_maneuver=maneuver_time, AttitudeRepresentation=SimDataRep);

% Three controllers:
%  (1) SD-OPT w/ shrinking horizon MPC at the terminal condition (used by default).
%  (2) SD-OPT w/ pure preview tracking - the typical online case.
%  (3) SD-MPC - full state-dependent MPC: full recursion over the costate and cost-to-go, at each timestep. The end of the window is treated via terminal cost Qyf, instead of the continued-to-be-held reference assumption ("SDOPT seed").

ctrl_fcn_1 = @(t,x,k,u_prev) compute_u_SDOPT(t, xmap(x), k, u_prev, ...
    r_, C, Qy, R, Qyf, qp, PreviewHorizon = SDOPTPreviewHorizon, DecimationOpts = dopts, ...
    SDC_A_function = A_fcn, ...
    SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

ctrl_fcn_2 = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, u_prev, ...
    r_, C, Qy, R, Qyf, qp, PreviewHorizon = SDOPTPreviewHorizon, DecimationOpts = dopts, ...
    AlwaysUseFullFiniteHorizonMPC = false, ...
    UseFullFiniteHorizonMPCAtTerminal = false, ...
    SDC_A_function = A_fcn, ...
    SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

ctrl_fcn_3 = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, u_prev, ...
    r_, C, Qy, R, Qyf, qp, PreviewHorizon = SDOPTPreviewHorizon, DecimationOpts = dopts, ...
    AlwaysUseFullFiniteHorizonMPC = true, ...
    SDC_A_function = A_fcn, ...
    SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

thrust_fcn = @(t) ones(1,6);
[~, ~, C_nat, Qy_nat, Qyf_nat] = get_weights(qp, SDCAttRep); % Evaluate cost in Native coordinates

ctrl_fcns = {ctrl_fcn_1, ctrl_fcn_2, ctrl_fcn_3};
labels = {'SD-OPT (MPC at terminal)', 'SD-OPT (pure preview)', 'SD-MPC'};

fprintf('\n--- SD-OPT decimated (with %s attitude) ---\n\n\n', SDCAttRep)
c1_     = cell(1,3);
t_      = cell(1,3);
X_      = cell(1,3);
U_      = cell(1,3);
stats_  = cell(1,3);
J_      = cell(1,3);
Jy_     = cell(1,3);
Ju_     = cell(1,3);
Jy_i_   = cell(1,3);
U_raw_  = cell(1,3);
for i=1:3
    [t, X, U, c1, c2, stats, U_raw] = run_sim(qp, tspan, x0, thrust_fcn, ctrl_fcns{i}, AttitudeRepresentation=SimDataRep);
    [J, Jy, Ju, Jy_i] = compute_J_LQT(t, X, U, Qy_nat, Qyf_nat, R, C_nat, ...
        x_ref_fcn, qp.Ts, AttitudeRepresentation=SimDataRep);
        
    % Store results
    t_{i} = t; X_{i} = X; U_{i} = U;
    c1_{i} = c1; stats_{i} = stats; U_raw_{i} = U_raw;
    J_{i} = J; Jy_{i} = Jy; Ju_{i} = Ju; Jy_i_{i} = Jy_i;

    fprintf('Config %d: %s \n', i, labels{i})
    fprintf('  Compute : %.2f s total → %.1f Hz equivalent\n', sum(c1), length(U)/sum(c1))
    fprintf('    running the simulation (rk4 or ode45) took up the other %.2f s\n', sum(c2))
    print_tracking_metrics(X, r_, J, Jy, Ju, Jy_i)
    print_saturation(U_raw, qp)
    fprintf('\n\n')
end


%% Compute and tracking comparison across the three configs

nr   = size(r_,2);
cmap = colororder;

figure(4); clf

ax_p = gobjects(1,3);
ax_u = gobjects(1,3);
ax_c = gobjects(1,3);
for i = 1:3

    t = t_{i}(1:nr);
    X = X_{i}(1:nr,:);
    U = U_{i};

    % --- Left: position vs reference ---

    ax_p(i) = subplot(3,3,3*i-2); hold on

    plot(t,X(:,1))
    plot(t,X(:,2))
    plot(t,X(:,3))

    plot(t,r_(1,:),'Color',cmap(1,:),'LineStyle','--')
    plot(t,r_(2,:),'Color',cmap(2,:),'LineStyle','--')
    plot(t,r_(3,:),'Color',cmap(3,:),'LineStyle','--')

    grid on
    yline(0,'Color',[1 1 1]*0.25)
    ylabel('Position [m]')
    legend('x','y','z')
    title(labels{i})

    % --- Middle: actuator inputs ---

    ax_u(i) = subplot(3,3,3*i-1); hold on

    % Actuator limits + nominal
    yline(qp.max_du,'Color',[1 1 1]*0.5)
    yline(0,'Color',[1 1 1]*0.67)
    yline(-qp.nominal_omegas,'Color',[1 1 1]*0.5)

    plot(t(1:end-1), U(1:nr-1,:)) % control input

    grid on
    ylabel('\delta U')

    % --- Right: compute breakdown ---

    t1 = [stats_{i}.TimeToComputeDiscreteSDCMatrices];
    t2 = t1 + [stats_{i}.TimeToSolveDARE];
    t3 = t2 + [stats_{i}.TimeToComputeFeedforward];

    dt1 = 1000*t1;
    dt2 = 1000*(t2-t1);
    dt3 = 1000*(t3-t2);

    ax_c(i) = subplot(3,3,3*i);
    area(t_{i}(1:end-1), [dt1(:) dt2(:) dt3(:)], 'EdgeColor', 'none')
    yline(qp.Ts * 1000, '--')
    ylabel("Compute [ms]")
    legend('sdc','dare','ff','Location','west')
end

links = [linkprop([ax_p ax_u ax_c], 'XLim'), ...
         linkprop(ax_p, 'YLim'), ...
         linkprop(ax_u, 'YLim'), ...
         linkprop(ax_c, 'YLim')];
setappdata(gcf, 'axis_links', links)  % links stay active only while referenced

xlim(ax_p(1), [0 maneuver_time])
ylim(ax_p(1), [-1 1]*1.1*max(abs(r_),[],'all'))
ylim(ax_u(1), [-1 1]*400)
ylim(ax_c(1), [-0.02 1]*qp.Ts*1.1 * 1000)

xlabel(ax_p(3), 'Time [s]')
xlabel(ax_u(3), 'Time [s]')
xlabel(ax_c(3), 'Time [s]')


%% Decimation comparison: SD-OPT (pure preview) and SD-MPC, with and without decimation

% Compares undecimated vs decimated: SD-OPT and SD-MPC. Print-only output.
% 
% NB: this section will take 1-2mins to run, as the undecimated MPC recursion at 1000Hz is very slow - expect it to take 1-2mins total

dec_labels = {'SD-OPT (pure preview), undecimated', 'SD-OPT (pure preview), decimated', 'SD-MPC, undecimated', 'SD-MPC, decimated'};
dec_opts   = {[], dopts, [], dopts};
full_mpc   = [false false true true];

N = length(r_);

fprintf('\n--- Decimation comparison (%g Hz, %.1fs preview, %s attitude) ---\n', 1/qp.Ts, SDOPTPreviewHorizon, SDCAttRep)
fprintf('  %-34s %8s %10s %8s %14s %10s %8s\n', 'Config', 'J', 'RMSE', '||u||', 'Total compute', 'Per-step', 'Budget')
for i = 1:4
    fprintf('  %-34s ', dec_labels{i});  % label first, so the user sees what's running
    ctrl_fcn = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, u_prev, ...
        r_, C, Qy, R, Qyf, qp, PreviewHorizon = SDOPTPreviewHorizon, DecimationOpts = dec_opts{i}, ...
        AlwaysUseFullFiniteHorizonMPC = full_mpc(i), ...
        UseFullFiniteHorizonMPCAtTerminal = false, ...
        SDC_A_function = A_fcn, ...
        SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));

    [t, X, U, c1] = run_sim(qp, tspan, x0, thrust_fcn, ctrl_fcn, AttitudeRepresentation=SimDataRep);
    J = compute_J_LQT(t, X, U, Qy_nat, Qyf_nat, R, C_nat, ...
        x_ref_fcn, qp.Ts, AttitudeRepresentation=SimDataRep);

    e_pos = X(1:N, 1:3) - r_(1:3, 1:N)';
    rmse  = 100 * sqrt(mean(e_pos.^2, 'all'));
    u_rms = sqrt(mean(sum(U(1:N-1,:).^2, 2)));  % RMS over time of the input vector norm
    ct    = c1(round(0.1/qp.Ts):end);  % discard noisy first samples

    fprintf('%8.2f %10s %8.1f %14s %10s %7.1f%%\n', J, sprintf('%.2f cm', rmse), u_rms, ...
        sprintf('%.2f s', sum(c1)), sprintf('%.3f ms', 1000*median(ct)), 100*median(ct)/qp.Ts)
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