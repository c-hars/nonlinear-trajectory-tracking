clear
root = fullfile(fileparts(mfilename('fullpath')), '..', '..'); addpath(root)
load_paths
load_copter_params

maxNumCompThreads(1);  % small matrices: multithreaded BLAS only adds timing jitter

% Loop-rate sweep: exact vs decimated, for SD-OPT (pure preview) and SD-MPC.
% 
%   Figure 1: per-step compute (median) vs loop rate, with the real-time threshold (1*Ts) visible.
%   Figure 2: tracking cost vs loop rate - the delta J of each decimated config (relative to exact's J).
% 
% The reference is padded past its end so that median compute remains meaningful (neither controller shrinks its horizon).
% 
% SD-OPT (MPC at terminal) is left out: its per-step compute sits between the SD-OPT and SD-MPC curves (preview control, then shrinking-horizon MPC at the terminal), and with the padding it would never reach its terminal phase anyway - it would just repeat SD-OPT.
% 
% Decimation uses one of two layouts: geometric or two-tier (see make_dopts).
% 
% NB: these are PC timings. On embedded hardware without optimised matrix routines and instruction sets, compute time follows the flop count much more closely; decimation there is where the benefit is clearest (and most useful).


%% Main parameters

SweepRates = [20 50 100 200 400 1000]; % [Hz]
maneuver_time = 5.0;       % [seconds]
SDOPTPreviewHorizon = 2.0; % [seconds]
SDCAttRep = 'Quaternion';  % Euler, Quaternion, MRP, FRA

PlotUnits = '%';
% 'ms': per-step compute [ms], against the sample period T_s
%  '%': per-step compute as a percentage of T_s


%% Loop-rate sweep

% Colour = controller, marker = decimation layout, dotted = held-x. Exact runs first per controller (the printed ΔJ needs its J).
cfgs = struct( ...
    'label',    {'SD-OPT, exact', 'SD-OPT, geometric', 'SD-OPT, two-tier', ...
                 'SD-MPC, exact', 'SD-MPC, geometric', 'SD-MPC, geometric w/ held-x', 'SD-MPC, two-tier', 'SD-MPC, two-tier w/ held-x'}, ...
    'full_mpc', {false, false, false, true, true, true, true, true}, ...
    'layout',   {'', 'geometric', 'two-tier', '', 'geometric', 'geometric', 'two-tier', 'two-tier'}, ...
    'held_x',   {false, false, false, false, false, true, false, true}, ...
    'style',    {'-o', '--s', '-.^', '-o', '--s', ':s', '-.^', ':^'});
nc       = numel(cfgs);
cidx     = 1 + [cfgs.full_mpc];
is_exact = cellfun(@isempty, {cfgs.layout});
i_ref    = arrayfun(@(c) find(is_exact & [cfgs.full_mpc] == c.full_mpc), cfgs);  % ΔJ baseline for each config: the undecimated run of the same controller
per_step = nan(numel(SweepRates), nc);
J        = nan(numel(SweepRates), nc);
cmap     = colororder;

fprintf('\n--- Loop-rate sweep (%.1fs maneuver, %.1fs preview, %s attitude) ---\n', maneuver_time, SDOPTPreviewHorizon, SDCAttRep)

for ir = 1:numel(SweepRates)
    rate  = SweepRates(ir);
    qp.Ts = 1/rate;
    prob  = setup_problem(qp, maneuver_time, SDCAttRep);
    N     = length(prob.r_);
    prob.r_ = [prob.r_, repmat(prob.r_(:,end), 1, round(SDOPTPreviewHorizon/qp.Ts))];  % hold the last reference value for one horizon past the end: SD-MPC keeps its full horizon throughout (as SD-OPT already does internally)

    fcns = cell(1, nc);
    for ic = 1:nc
        if is_exact(ic), dopts = []; else, dopts = make_dopts(qp.Ts, cfgs(ic).layout, cfgs(ic).held_x); end
        fcns{ic} = make_ctrl(prob, SDOPTPreviewHorizon, cfgs(ic).full_mpc, dopts);
    end

    % Untimed warm-up: mitigates JIT and CPU clock ramp-up skewing the timings
    for ic = 1:nc
        for k = 1:min(100, N-1)
            fcns{ic}((k-1)*qp.Ts, prob.x0, k, zeros(6,1));
        end
    end

    % Each entry is printed as soon as its run finishes
    fprintf('  %g Hz\n', rate)
    for ic = 1:nc
        [t, X, U, ~, ~, stats] = run_sim(qp, prob.tspan, prob.x0, prob.thrust_fcn, fcns{ic}, ...
            AttitudeRepresentation=prob.SimDataRep);
        per_step(ir,ic) = compute_per_step(stats, qp.Ts);
        J(ir,ic) = compute_J_LQT(t, X, U, prob.Qy, prob.Qyf, prob.R, prob.C, prob.x_ref_fcn, qp.Ts, ...
            AttitudeRepresentation=prob.SimDataRep);
        fprintf('    %-28s %9s  %7s', cfgs(ic).label, sprintf('%.3f ms', 1000*per_step(ir,ic)), ...
            sprintf('(%.1f%%)', 100*per_step(ir,ic)/qp.Ts))
        if ~is_exact(ic)
            fprintf('   ΔJ %+.2f%%', 100*(J(ir,ic)/J(ir,i_ref(ic)) - 1))
        end
        fprintf('\n')
    end

    % Update the plot as results come in
    rates = SweepRates(1:ir);
    switch PlotUnits
        case 'ms'
            y       = 1000 * per_step(1:ir,:);
            y_limit = 1000 ./ rates;
            y_label = 'Per-step compute [ms, median]';
        case '%'
            y       = 100 * per_step(1:ir,:) .* rates(:);
            y_limit = 100 * ones(size(rates));
            y_label = 'Per-step compute [% of T_s, median]';
        otherwise
            error("Unknown PlotUnits: '%s'", PlotUnits)
    end

    dJ = 100 * (J(1:ir,:) ./ J(1:ir,i_ref) - 1);  % [%] vs exact, same controller

    % Figure 1: per-step compute
    figure(1); clf; hold on
    h = gobjects(1, nc+1);
    for ic = 1:nc
        h(ic) = plot(rates, y(:,ic), cfgs(ic).style, 'Color', cmap(cidx(ic),:));
    end
    h(nc+1) = plot(rates, y_limit, '--', 'Color', [1 1 1]*0.25);
    set(gca, 'XScale', 'log', 'YScale', 'log')
    grid on
    xlim([SweepRates(1) SweepRates(end)])
    ylim([0.8*min(y,[],'all') 1.5*max(y,[],'all')])  % fit the data; the real-time limit line can run off the top
    xlabel('Loop rate [Hz]')
    ylabel(y_label)
    leg   = [{cfgs.label}, {'real-time limit (T_s)'}];
    order = [1 4 9, 2 5 6, 3 7 8];  % legend columns: exact | geometric | two-tier; rows: SD-OPT, SD-MPC, held-x (indices into cfgs)
    legend(h(order), leg(order), 'Location', 'northwest', 'NumColumns', 3)
    title({'Loop rate vs per-step compute', ...
            sprintf('\\fontsize{10}%g–%gHz | %.1f s prediction horizon', SweepRates(1), SweepRates(end), SDOPTPreviewHorizon)})

    % Figure 2: tracking cost, relative to exact
    figure(2); clf; hold on
    h = gobjects(1, nc);
    for ic = find(~is_exact)
        h(ic) = plot(rates, dJ(:,ic), cfgs(ic).style, 'Color', cmap(cidx(ic),:));
    end
    yline(0, 'Color', [1 1 1]*0.25)
    set(gca, 'XScale', 'log')
    grid on
    xlim([SweepRates(1) SweepRates(end)])
    xlabel('Loop rate [Hz]')
    ylabel('\DeltaJ vs exact [%]')
    order = [2 5 6, 3 7 8];  % as Figure 1, without the exact column
    legend(h(order), {cfgs(order).label}, 'Location', 'northwest', 'NumColumns', 2)
    title({'Loop rate vs tracking cost', ...
         sprintf('\\fontsize{10}%g–%gHz | %.1f s prediction horizon', SweepRates(1), SweepRates(end), SDOPTPreviewHorizon)})
    drawnow
end


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

function dopts = make_dopts(Ts, layout, held_x)
% Decimation layouts, both reaching a ~20Hz far-horizon preview:
%   geometric: factors 1, 2, 4, ... up to ~20Hz
%   two-tier:  factor 1, then straight to ~20Hz
% Geometric adds a near tier per doubling of the rate, so its recursion grows faster than two-tier's: at 1000 Hz (2 s preview), 206 vs 139 steps (undecimated: 1999).
    n = floor(log2((1/20)/Ts));  % far-horizon factor: 2^n
    switch layout
        case 'geometric'
            df = 2.^(0:n);
            dopts = struct('DecimationFactors', df, 'StepsPerTier', df(end), 'UseHeldXAssumption', held_x);
        case 'two-tier'
            df = unique(2.^[0, n]);  % unique: just [1] at 20Hz
            dopts = struct('DecimationFactors', df, 'StepsPerTier', 2*df(end), 'UseHeldXAssumption', held_x);
        otherwise
            error("Unknown layout: '%s'", layout)
    end
end

function fcn = make_ctrl(prob, horizon, full_mpc, dopts)
% SD-OPT (pure preview): full_mpc = false
% SD-MPC:                full_mpc = true
% dopts:  decimation options ([] = exact)
    xmap = prob.xmap;
    fcn = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, xmap(x), k, u_prev, ...
        prob.r_, prob.C, prob.Qy, prob.R, prob.Qyf, prob.qp, ...
        PreviewHorizon = horizon, ...
        DecimationOpts = dopts, ...
        AlwaysUseFullFiniteHorizonMPC = full_mpc, ...
        UseFullFiniteHorizonMPCAtTerminal = false, ...
        SDC_A_function = prob.A_fcn, ...
        SDC_B_function = @(uk,qp) get_B_matrix_SDRE(uk,qp));
end

function per_step = compute_per_step(stats, Ts)
% Median controller compute per step, from the solve_info breakdown
    ct = [stats.TimeToComputeDiscreteSDCMatrices] + [stats.TimeToSolveDARE] + [stats.TimeToComputeFeedforward];
    per_step = median(ct(round(0.1/Ts):end));  % discard noisy first samples
end
