clear
load_paths
load_copter_params

maxNumCompThreads(1);  % small matrices: multithreaded BLAS only adds timing jitter

% Decimated recursions: exact (undecimated) vs decimated, for one controller config (chosen below).
%
%   Figure 1: the preview grid - how the preview horizon is laid out into tiers of fine and coarse steps.
%   Figure 2: closed-loop runs, one row per variant (exact, decimated, decimated w/ held-x): position, attitude, inputs, compute.
%   Figure 3: closed-loop differences vs exact (decimated - exact): position, attitude, inputs, each variant on its own trajectory.
%   Figure 4: exact vs decimated in closed loop: position error and per-step compute, each variant on its own trajectory.
%   Figure 5: replay along the exact trajectory: the difference in control input when each decimated controller is evaluated at the exact controller's states.
%
% Decimation applies to the preview recursion (SD-OPT) and to the full Riccati recursion (SD-MPC). See docs/decimation/decimation.md for the details.
% Per-step compute vs loop rate is covered separately, in demo_decimation_rates.m.

% NB: these are PC timings. On embedded hardware without optimised matrix routines and instruction sets, compute time follows the flop count much more closely; decimation there is where the benefit is clearest (and most useful).


%% Controller config: choose one
%
%  'SD-OPT (MPC at terminal)': preview control, switching to shrinking-horizon SD-MPC once the end of the reference comes within the preview horizon. This is the default SD-OPT setup (as in compute_u_SDOPT).
%  'SD-OPT (pure preview)':    preview control throughout, with the reference held past its end - the typical online case (no shrinking-horizon effects at the terminal).
%  'SD-MPC':                   full state-dependent MPC: the full Riccati recursion over the cost-to-go and costate, every step. The end of the window is treated via terminal cost Qyf, instead of the held-reference assumption ("SDOPT seed").
%
% Controllers with an MPC recursion (SD-MPC, and the terminal phase of SD-OPT (MPC at terminal)) also get held-x on/off variants - see UseHeldXAssumption, DecimationOpts, compute_u_SDDRE_v3.

Config = 'SD-OPT (MPC at terminal)';
% Config = 'SD-OPT (pure preview)';
% Config = 'SD-MPC';


%% Main parameters

qp.Ts = 1/1000;            % sample rate
maneuver_time = 5.0;       % [seconds] - 5s is aggressive but non-saturating: differences reflect the control laws rather than clipping
SDOPTPreviewHorizon = 2.0; % [seconds]
SDCAttRep = 'Quaternion';  % Euler, Quaternion, MRP, FRA
InputUnits = 'RPM';        % RPM, rad/s - inputs are plotted as the deviation from hover either way

prob  = setup_problem(qp, maneuver_time, SDCAttRep);
dopts = make_dopts(qp.Ts);
dopts_hx = dopts; dopts_hx.UseHeldXAssumption = true;

switch Config
    case 'SD-OPT (MPC at terminal)', full_mpc = false; mpc_at_terminal = true;
    case 'SD-OPT (pure preview)',    full_mpc = false; mpc_at_terminal = false;
    case 'SD-MPC',                   full_mpc = true;  mpc_at_terminal = false;
    otherwise, error("Unknown Config: '%s'", Config)
end

% Variants, in run order.
% Decimated first, exact last (it's the slow one). Exact also evaluates the decimated controllers at its own states (for Figure 5).
% (cidx: colour index, consistent across the figures)
if full_mpc || mpc_at_terminal
    variants = struct('label',  {'decimated (w/ held-x)', 'decimated', 'exact'}, ...
                      'dopts',  {dopts_hx, dopts, []}, ...
                      'cidx',   {4, 2, 1});
else
    variants = struct('label',  {'decimated', 'exact'}, ...
                      'dopts',  {dopts, []}, ...
                      'cidx',   {2, 1});
end
nv = numel(variants);

N         = length(prob.r_);
n_preview = round(SDOPTPreviewHorizon/qp.Ts);
k_switch  = max(N - n_preview, 1);   % SD-OPT (MPC at terminal): first step of the terminal MPC phase
t_switch  = (k_switch - 1) * qp.Ts;

w_despike = 9;  % moving-min window [samples] for compute: removes OS spikes (tuned at 1000 Hz)

switch InputUnits  % U is in rad/s
    case 'RPM',   u_scale = 60/(2*pi);
    case 'rad/s', u_scale = 1;
    otherwise, error("Unknown InputUnits: '%s'", InputUnits)
end


%% Figure 1: the preview grid

% One controller call at the initial state, just to get the grid it builds
f_dec = make_ctrl(prob, SDOPTPreviewHorizon, full_mpc, mpc_at_terminal, dopts);
[~, info0] = f_dec(0, prob.x0, 1, zeros(6,1));
grid_ = info0.PreviewGrid;
n_rec = grid_.n_fine_rem + sum(grid_.n_blocks);

fprintf('\n--- Preview grid (%g Hz, %.1fs preview) ---\n', 1/qp.Ts, SDOPTPreviewHorizon)
fprintf('  Decimation factors: %s, %d blocks per near tier\n', mat2str(grid_.df_), grid_.n_steps_per_tier)
fprintf('  %d fine steps → %d recursion steps (%.1fx reduction)\n', grid_.n_fine, n_rec, grid_.n_fine/n_rec)

figure(1); clf

ax1 = subplot(2,1,1);
h = draw_grid_strip(grid_, qp.Ts);
set(ax1, 'XScale', 'log')
xlim([qp.Ts SDOPTPreviewHorizon])
leg = [arrayfun(@(d) sprintf('d = %d', d), grid_.df_, 'UniformOutput', false), {'fine remainder'}];
valid = isgraphics(h);
legend(h(valid), leg(valid), 'Location', 'northoutside', 'NumColumns', min(nnz(valid), 4))

ax2 = subplot(2,1,2);
draw_grid_strip(grid_, qp.Ts);
xlim([0 SDOPTPreviewHorizon])
xlabel('Preview time [s]')

sgtitle({sprintf('Preview grid: %d fine steps → %d recursion steps', grid_.n_fine, n_rec), ...
         sprintf('\\fontsize{10}%g Hz, %.1f s preview', 1/qp.Ts, SDOPTPreviewHorizon)})
drawnow


%% Closed-loop runs: decimated, then exact

fcns = arrayfun(@(v) make_ctrl(prob, SDOPTPreviewHorizon, full_mpc, mpc_at_terminal, v.dopts), ...
    variants, 'UniformOutput', false);

fprintf('\n--- %s (%.1fs maneuver, %g Hz, %.1fs preview, %s attitude) ---\n', ...
    Config, prob.tspan(end), 1/qp.Ts, SDOPTPreviewHorizon, SDCAttRep)

% SD-OPT (MPC at terminal): per-step median over the preview phase only.
% In the terminal phase, per-step cost isn't steady - it falls steadily as the horizon shrinks.
% The two decimated variants run identical code until the terminal switch, so their medians should match: any gap is measurement noise, a rough noise floor for the timing comparisons.
if mpc_at_terminal
    k_med_end = k_switch - 1;
    med_label = 'Per-step* (med)';
else
    k_med_end = Inf;
    med_label = 'Per-step (med)';
end

fprintf('  %-24s  %8s  %10s  %8s  %13s  %15s  %14s\n', 'Variant', 'J', 'RMSE', '||u||', 'Total compute', ...
    med_label, 'Per-step (max)')

% Untimed warm-up: mitigates JIT and CPU clock ramp-up skewing the timings
k_warm = [1:100, k_switch:min(k_switch+99, N-1)];
for i = 1:nv
    for k = k_warm
        fcns{i}((k-1)*qp.Ts, prob.x0, k, zeros(6,1));
    end
end

% One closed-loop run per variant
res = cell(1, nv);
for i = 1:nv
    if i < nv
        % Decimated: run as is
        ctrl = fcns{i};
    else
        % Exact (last): the decimated controllers shadow it - evaluated at each of its states, outputs logged (Figure 5) but not applied
        ctrl = @(t,x,k,u_prev) with_shadows(fcns{nv}, fcns(1:nv-1), t, x, k, u_prev);
    end

    fprintf('  %-24s   0%%', variants(i).label)
    ctrl_p = @(t,x,k,u_prev) with_progress(ctrl, t, x, k, u_prev, N-1);
    [t, X, U, ~, ~, stats] = run_sim(qp, prob.tspan, prob.x0, prob.thrust_fcn, ctrl_p, AttitudeRepresentation=prob.SimDataRep);
    fprintf('\b\b\b\b')

    res{i} = metrics(t, X, U, stats, prob, w_despike, k_med_end);
    print_row(res{i})
end
if mpc_at_terminal
    fprintf('  (*) Over the preview phase only (median only makes sense here).\n')
end


%% Figure 2: each variant's closed loop run - position, attitude, inputs, compute 

rows = nv:-1:1;  % exact on top
cmap = colororder;

figure(2); clf

ax_p = gobjects(1,nv);
ax_a = gobjects(1,nv);
ax_u = gobjects(1,nv);
ax_c = gobjects(1,nv);
att_max = 0;
for r = 1:nv
    i = rows(r);
    t = res{i}.t(1:N);
    X = res{i}.X(1:N,:);
    U = res{i}.U;

    % --- Position vs reference ---

    ax_p(r) = subplot(nv,4,4*r-3); hold on

    plot(t,X(:,1))
    plot(t,X(:,2))
    plot(t,X(:,3))

    plot(t,prob.r_(1,:),'Color',cmap(1,:),'LineStyle','--')
    plot(t,prob.r_(2,:),'Color',cmap(2,:),'LineStyle','--')
    plot(t,prob.r_(3,:),'Color',cmap(3,:),'LineStyle','--')

    grid on
    yline(0,'Color',[1 1 1]*0.25)
    ylabel('Position [m]')
    title(variants(i).label)
    if r == 1, legend('x','y','z'), end

    % --- Attitude ---

    ax_a(r) = subplot(nv,4,4*r-2); hold on

    Xe = to_euler_state(X, prob.SimDataRep);
    plot(t, Xe(:,7:9)*180/pi)
    att_max = max(att_max, max(abs(Xe(:,7:9)*180/pi), [], 'all'));

    grid on
    yline(0,'Color',[1 1 1]*0.25,'HandleVisibility','off')
    ylabel('Euler angles [deg]')
    if r == 1, legend('\phi','\theta','\psi'), end

    % --- Actuator inputs ---

    ax_u(r) = subplot(nv,4,4*r-1); hold on

    % Actuator limits + nominal
    yline(u_scale*qp.max_du,'Color',[1 1 1]*0.5)
    yline(0,'Color',[1 1 1]*0.67)
    yline(-u_scale*qp.nominal_omegas,'Color',[1 1 1]*0.5)

    plot(t(1:end-1), u_scale*U(1:N-1,:)) % control input

    grid on
    ylabel(sprintf('U [%s]', InputUnits))

    % --- Compute breakdown ---

    dt1 = movmin(1000*res{i}.ct_sdc,  w_despike);
    dt2 = movmin(1000*res{i}.ct_dare, w_despike);
    dt3 = movmin(1000*res{i}.ct_ff,   w_despike);

    ax_c(r) = subplot(nv,4,4*r);
    area(res{i}.t(1:numel(dt1)), [dt1(:) dt2(:) dt3(:)], 'EdgeColor', 'none')
    yline(qp.Ts * 1000, '--')
    if mpc_at_terminal
        xline(t_switch, ':', 'terminal MPC', 'LabelVerticalAlignment', 'top')
    end
    ylabel('Compute [ms]')
    if r == 1, legend('sdc','dare','ff / recursion','Location','west'), end
end

links = [linkprop([ax_p ax_a ax_u ax_c], 'XLim'), ...
         linkprop(ax_p, 'YLim'), ...
         linkprop(ax_a, 'YLim'), ...
         linkprop(ax_u, 'YLim'), ...
         linkprop(ax_c, 'YLim')];
setappdata(gcf, 'axis_links', links)  % keep the link objects alive with the figure (linkprop links stop working once their object is deleted, e.g. by clear)

% Compute axis: fit all variants, but cap at 2x the budget (exact MPC can run far off the top)
yl_upper = min(2*qp.Ts, 1.25*max(cellfun(@(m) prctile(m.ct, 95), res))) * 1000;

xlim(ax_p(1), [0 prob.tspan(end)])
ylim(ax_p(1), [-1 1]*1.1*max(abs(prob.r_),[],'all'))
ylim(ax_a(1), [-1 1]*1.1*att_max)
ylim(ax_u(1), [-1 1]*400*u_scale)
ylim(ax_c(1), [-0.01 1]*yl_upper)

xlabel(ax_p(nv), 'Time [s]')
xlabel(ax_a(nv), 'Time [s]')
xlabel(ax_u(nv), 'Time [s]')
xlabel(ax_c(nv), 'Time [s]')

sgtitle({'Decimated vs exact: states, inputs and compute', ...
         sprintf('\\fontsize{10}%s, %g Hz, %.1f s preview', Config, 1/qp.Ts, SDOPTPreviewHorizon)})
drawnow


%% Figure 3: closed-loop differences vs exact - position, attitude, inputs

figure(3); clf
Xe_ex = to_euler_state(res{nv}.X(1:N,:), prob.SimDataRep);
ax_p = gobjects(1,nv-1);
ax_a = gobjects(1,nv-1);
ax_u = gobjects(1,nv-1);
for r = 1:nv-1
    i = nv - r;  % decimated, then decimated (w/ held-x)
    t = res{i}.t(1:N);
    dX = res{i}.X(1:N,:) - res{nv}.X(1:N,:);
    dU = u_scale * (res{i}.U(1:N-1,:) - res{nv}.U(1:N-1,:));

    % --- Position ---

    ax_p(r) = subplot(nv-1,3,3*r-2); hold on
    plot(t, 100*dX(:,1:3))
    grid on
    ylabel('\Deltaposition [cm]')
    title(variants(i).label)
    if r == 1, legend('x','y','z'), end

    % --- Attitude ---

    ax_a(r) = subplot(nv-1,3,3*r-1); hold on
    Xe = to_euler_state(res{i}.X(1:N,:), prob.SimDataRep);
    dE = mod((Xe(:,7:9) - Xe_ex(:,7:9))*180/pi + 180, 360) - 180;  % wrap: a heading wraparound isn't a 360 deg difference
    plot(t, dE)
    grid on
    ylabel('\DeltaEuler angles [deg]')
    if r == 1, legend('\phi','\theta','\psi'), end

    % --- Actuator inputs ---

    ax_u(r) = subplot(nv-1,3,3*r); hold on
    plot(t(1:end-1), dU)
    grid on
    ylabel(sprintf('\\DeltaU [%s]', InputUnits))
end

linkaxes([ax_p ax_a ax_u], 'x')
xlim(ax_p(1), [0 prob.tspan(end)])
xlabel(ax_p(nv-1), 'Time [s]')
xlabel(ax_a(nv-1), 'Time [s]')
xlabel(ax_u(nv-1), 'Time [s]')

sgtitle({'Decimated vs exact: state and input differences', ...
         sprintf('\\fontsize{10}\\Delta = decimated - exact · %s, %g Hz, %.1f s preview', Config, 1/qp.Ts, SDOPTPreviewHorizon)})
drawnow


%% Figure 4: exact vs decimated in closed loop - position error, per-step compute

figure(4); clf
ax = gobjects(1,2);

% --- Position tracking error ---

ax(1) = subplot(2,1,1); hold on
for i = rows
    ls = '-';
    if i < nv, ls = '--'; end
    plot(res{i}.t(1:N), res{i}.e_pos, 'Color', cmap(variants(i).cidx,:), 'LineStyle', ls)
end
grid on
ylabel('Position error [cm]')
legend({variants(rows).label})

% --- Per-step compute ---

ax(2) = subplot(2,1,2); hold on
for i = rows
    plot(res{i}.t(1:numel(res{i}.ct)), 1000*res{i}.ct, 'Color', cmap(variants(i).cidx,:))
end
yline(qp.Ts * 1000, '--')
if mpc_at_terminal
    xline(t_switch, ':', 'terminal MPC', 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off')
end
grid on
ylabel('Compute [ms]')
legend([{variants(rows).label}, {'T_s'}], 'Location', 'northwest')
ylim([0 1.25*max(cellfun(@(m) prctile(m.ct, 95), res))*1000])

linkaxes(ax, 'x')
xlim(ax(1), [0 prob.tspan(end)])
xlabel('Time [s]')

sgtitle({'Decimated vs exact: tracking error and compute', ...
         sprintf('\\fontsize{10}%s, %g Hz, %.1f s preview', Config, 1/qp.Ts, SDOPTPreviewHorizon)})
drawnow


%% Figure 5: control differences at the same state (replay along exact's trajectory)
% Each decimated controller sees exact's inputs (x_k, u_{k-1}) at every step, so any difference is down to the control law alone.

Um = vertcat(res{nv}.stats.UMain);
Us = cat(3, res{nv}.stats.UShadow);  % (decimated variant) x (motor) x (step)
for i = 1:nv-1
    dU = u_scale * (reshape(Us(i,:,:), size(Us,2), [])' - Um);  % [InputUnits], step x motor
    res{i}.du_2   = vecnorm(dU, 2, 2);
    res{i}.du_inf = max(abs(dU), [], 2);
end

figure(5); clf
ax = gobjects(1,nv-1);
for r = 1:nv-1
    i = nv - r;  % decimated, then decimated (w/ held-x)
    col = cmap(variants(i).cidx,:);
    t_r = res{nv}.t(1:numel(res{i}.du_2));

    ax(r) = subplot(nv-1,1,r); hold on
    plot(t_r, res{i}.du_2,   'Color', col)
    plot(t_r, res{i}.du_inf, 'Color', 0.6*col, 'LineStyle', '--')
    if mpc_at_terminal
        xline(t_switch, ':', 'terminal MPC', 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off')
    end
    grid on
    ylabel(sprintf('\\Deltau [%s]', InputUnits))
    legend('2-norm (all motors)', '\infty-norm (worst motor)', 'Location', 'northwest')
    title(variants(i).label)
end

linkaxes(ax, 'x')
xlim(ax(1), [0 prob.tspan(end)])
xlabel('Time [s]')

sgtitle({'Decimated vs exact: control discrepancy \Deltau', ...
         sprintf('\\fontsize{10}each an open-loop replay against exact''s \\{x_k, u_{k-1}\\} trajectory · %s, %g Hz, %.1f s preview', Config, 1/qp.Ts, SDOPTPreviewHorizon)})
drawnow


%% ---- Local functions ----

function prob = setup_problem(qp, maneuver_time, SDCAttRep)
% Weights, reference, and coordinate map for a given sample rate.
% get_weights() rescales Qy, Qyf with Ts, so this is rebuilt per rate.

    prob.qp         = qp;
    prob.SimDataRep = 'Quaternion';
    prob.xmap       = xq2sdc(SDCAttRep);
    prob.A_fcn      = get_SDC_A_function(SDCAttRep);
    prob.thrust_fcn = @(t) ones(1,6);

    [~, prob.R, prob.C, prob.Qy, prob.Qyf] = get_weights(qp, SDCAttRep);
    [prob.r_, prob.x_ref_fcn, ~, prob.tspan, prob.x0] = load_fig8_traj(qp, prob.C, ...
        t_maneuver=maneuver_time, AttitudeRepresentation=prob.SimDataRep);
end

function dopts = make_dopts(Ts)
% Geometric progression of decimation factors (1, 2, 4, ...), up to a ~20Hz far-horizon preview.
% This is the default layout, but the decimation options are an open design parameter.
    df = 2.^(0:floor(log2((1/20)/Ts)));
    dopts = struct('DecimationFactors', df, 'StepsPerTier', df(end), 'UseHeldXAssumption', false);

    % Alternative: two tiers - fine, then straight to ~20Hz (no intermediate geometric progression)
    % df = 2.^[0, floor(log2((1/20)/Ts))];
    % dopts = struct('DecimationFactors', df, 'StepsPerTier', 2*df(end), 'UseHeldXAssumption', false);
end

function fcn = make_ctrl(prob, horizon, full_mpc, mpc_at_terminal, dopts)
% SD-OPT (MPC at terminal): mpc_at_terminal
% SD-OPT (pure preview):    neither flag
% SD-MPC:                   full_mpc
% dopts:  decimation options ([] = exact)
    xmap = prob.xmap;
    fcn = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, u_prev, ...
        prob.r_, prob.C, prob.Qy, prob.R, prob.Qyf, prob.qp, ...
        PreviewHorizon = horizon, ...
        DecimationOpts = dopts, ...
        AlwaysUseFullFiniteHorizonMPC = full_mpc, ...
        UseFullFiniteHorizonMPCAtTerminal = mpc_at_terminal, ...
        SDC_A_function = prob.A_fcn, ...
        SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));
end

function [u, info] = with_shadows(main_fcn, shadow_fcns, t, x, k, u_prev)
% Run main_fcn in the loop, and also evaluate each of shadow_fcns at the same state for comparison
% Outputs are logged in info (UMain, UShadow); only main_fcn's u is applied.
    [u, info] = main_fcn(t, x, k, u_prev);
    info.UMain   = reshape(u, 1, []);
    info.UShadow = zeros(numel(shadow_fcns), numel(u));
    for j = 1:numel(shadow_fcns)
        info.UShadow(j,:) = shadow_fcns{j}(t, x, k, u_prev);
    end
end

function varargout = with_progress(ctrl_fcn, t, x, k, u_prev, Nsteps)
    if mod(k, ceil(Nsteps/10)) == 0
        fprintf('\b\b\b\b%3d%%', min(100, round(100*k/Nsteps)));
    end
    [varargout{1:nargout}] = ctrl_fcn(t, x, k, u_prev);
end

function m = metrics(t, X, U, stats, prob, w_despike, k_med_end)
    qp = prob.qp;
    N  = length(prob.r_);

    m.J = compute_J_LQT(t, X, U, prob.Qy, prob.Qyf, prob.R, prob.C, ...
        prob.x_ref_fcn, qp.Ts, AttitudeRepresentation=prob.SimDataRep);

    e_pos   = X(1:N, 1:3) - prob.r_(1:3, 1:N)';
    m.rmse  = 100 * sqrt(mean(e_pos.^2, 'all'));
    m.u_rms = sqrt(mean(sum(U(1:N-1,:).^2, 2)));  % RMS over time of the input vector norm
    m.e_pos = 100 * vecnorm(e_pos, 2, 2);

    % Controller compute per step, from the solve_info breakdown (excludes any shadow evaluations)
    k_settle  = round(0.1/qp.Ts);  % discard noisy first samples
    m.ct_sdc  = [stats.TimeToComputeDiscreteSDCMatrices];
    m.ct_dare = [stats.TimeToSolveDARE];
    m.ct_ff   = [stats.TimeToComputeFeedforward];
    m.ct      = m.ct_sdc + m.ct_dare + m.ct_ff;
    m.total    = sum(m.ct);
    m.per_step = median(m.ct(k_settle:min(k_med_end, end)));
    m.max_step = max(movmin(m.ct(k_settle:end), w_despike));  % despiked max

    m.t = t; m.X = X; m.U = U; m.stats = stats;
end

function print_row(m)
    fprintf(' %8.2f  %10s  %8.1f  %13s  %15s  %14s\n', m.J, sprintf('%.2f cm', m.rmse), m.u_rms, ...
        sprintf('%.2f s', m.total), sprintf('%.3f ms', 1000*m.per_step), sprintf('%.3f ms', 1000*m.max_step))
end

function h = draw_grid_strip(grid_, Ts)
% Preview grid as a strip, coloured by tier

    hold on
    cmap = colororder;
    h = gobjects(1, grid_.n_tiers + 1);
    grey = [1 1 1]*0.6;
    lw = 0.25;

    % Fine remainder
    for j = 1:grid_.n_fine_rem
        h(end) = patch(Ts*[j j+1 j+1 j], [0 0 1 1], grey, 'EdgeColor', 0.75*grey, 'LineWidth', lw);
    end
    tiers_x = zeros(0,2);  % tier extents, outlined last so no block edge is drawn over them
    if grid_.n_fine_rem > 0
        tiers_x(end+1,:) = Ts*[1, 1+grid_.n_fine_rem];
    end

    % Tiers
    for i = 1:grid_.n_tiers
        if grid_.n_blocks(i) == 0, continue, end
        d   = grid_.df_(i);
        col = cmap(mod(i-1, size(cmap,1)) + 1, :);
        for b = 1:grid_.n_blocks(i)
            j0 = grid_.tier_start(i) + (b-1)*d;
            h(i) = patch(Ts*[j0 j0+d j0+d j0], [0 0 1 1], col, 'EdgeColor', 0.75*col, 'LineWidth', lw);
        end
        tiers_x(end+1,:) = Ts*[grid_.tier_start(i), grid_.tier_start(i) + grid_.n_blocks(i)*d];
    end

    for r = 1:size(tiers_x,1)
        x = tiers_x(r,:);
        patch(x([1 2 2 1]), [0 0 1 1], 'k', 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 3*lw, 'HandleVisibility', 'off')
    end
    set(gca, 'SortMethod', 'childorder')  % draw in creation order (outlines on top)

    set(gca, 'YTick', [])
    grid on
    ylim([-0.2 1.2])
end
