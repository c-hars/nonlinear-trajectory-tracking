# Nonlinear near-optimal tracking: SDDRE control

**Nonlinear near-optimal tracking** with computational cost similar to linear MPC — a receding-horizon State-Dependent Riccati controller with reference preview, warm-started Riccati solves, and graceful degradation by construction.

Design philosophy: similar to nonlinear MPC, but with a predictable and bounded WCET and fixed memory requirements; the controller should extend linear tracking to handle nonlinear dynamics, but with a design that prioritises real-time suitability and embedded compute.
<table>

<!-- <tr>
  <th width="50%" align="center" style="text-align:center">Linear control (Linear Quadratic Tracking)</th>
  <th width="30%" align="center" style="text-align:center">SDDRE-based control</th>
</tr> -->

<tr>
  <td align="center">
    <img src="docs/figures/animation_ghost_4p5_enddelay.gif" width="636"><br>
    Demonstration: trajectory tracking in 3D space. The figure shows the trajectory under the LQT-based controller (dark grey, only visible when t≤1.5s) and the SDDRE-based controller (coloured; 3D orientation of the hexacopter[^1] indicated by the red/blue/yellow axes). The LQT-based controller quickly yields instability on this aggressive trajectory, with the feedback unable to fully mitigate the disturbance resulting from model mismatch (small-angle assumptions); performance loss and instability result. In contrast, the SDDRE-based controller successfully completes the maneuver; controls are adjusted via deterministic and efficiently-computed laws in response to the nonlinearities — effectively extending the LQT method to the nonlinear case.
  </td>
  <!-- <td align="center">
    <img src="docs/figures/figure8_linearAndSddre_4p5.png" width="500"><br>
  </td> -->
</tr>

</table>

[^1]: NB: the hexacopter is shown at 1.5× scale for improved visibility, which can make the trajectory appear smaller than it is; refer to the x/y/z axes for the true scale. Units: [m].

The implemented algorithm can be thought of as LQT-based MPC using SDC matrices, with additional design choices focused on improving robustness, tracking reliability, and efficient (NLP and optimiser-free) compute.

<table align="center">
  <tr>
    <th width="50%" align="center" style="text-align:center">Linear control (Linear Quadratic Tracking)</th>
    <th width="50%" align="center" style="text-align:center">SDDRE-based control</th>
  </tr>
  <tr>
      <td align="center">
      <img src="docs/figures/figure8_linear_8p5.png" width="475"><br>
      </td>
      <td align="center">
      <img src="docs/figures/figure8_sddre_8p5.png" width="475"><br>
      </td>
  </tr>
  <tr>
      <td align="center" colspan="2">
      <img src="docs/figures/errors_and_efforts.png" width="600"><br>
      Compared to the baseline linear method (Linear Quadratic Tracking, LQT), the SDDRE-based method improves tracking error, extends the usable range of actuator authority, and enables more aggressive maneuvering.
      </td>
  </tr>
</table>

---

## Summary

This repository implements a discrete-time **SDDRE tracking controller** for a 6DOF hexacopter, with benchmarks against finite-horizon Linear Quadratic Tracking (LQT) on a family of aggressive 3D maneuvers. At each control step the controller:

1. re-factorises the nonlinear dynamics into state-dependent coefficient (SDC) form $\dot x = A(x)\,x + B(u)\,u$ at the current state,
2. solves the resulting discrete algebraic Riccati equation, yielding the **infinite-horizon** solution $P_{ss}$ and gain $K_{ss}$ used for feedback — the solve for $P_{ss}$ is warm-started from the previous step's solution via Newton–Kleinman iteration (typically converging in 3-4 cheap Lyapunov solves),
3. computes a **preview feedforward** by a fixed-length backward costate recursion over an upcoming reference window, seeded so that the controller degrades *smoothly* to the infinite-horizon solution as preview information runs out,
4. switches to the exact finite-horizon Riccati recursion near the trajectory terminus, so the terminal cost is honoured exactly.

On a fixed-shape maneuver swept over duration $t_f$, the linear controller loses the trajectory below $t_f \approx 9$ s, while SDDRE tracks down to $t_f \approx 4.5$ s — roughly a 2× expansion of the feasible maneuver envelope, i.e. **capable of handling a ~4× higher acceleration demand**, along with general improvements in tracking error and control effort; and crucially, with per-step compute that is a *fixed* sequence of dense linear-algebra operations: no online NLP, no line searches, no variable-iteration solvers – consistent and bounded per-step compute.

---

## Motivation

Optimal tracking controllers face a familiar tradeoff:

- **Linear (LQR/LQT) designs** are cheap, certifiable, and well understood — but they are coupled to the hover linearisation. Once a maneuver demands sustained large attitude excursions, the small-angle model is simply wrong, and tracking collapses no matter how the weights are tuned.
- **Nonlinear MPC** handles the full dynamics and constraints — at the price of an online NLP whose solve time varies with iteration count, active-set changes, and initialisation quality. For a flight controller, *worst-case* timing is what matters, and it is hard to bound without RTI-class machinery and the engineering it entails. While backed by strong theory, practically NLMPC performance can be sensitive to prediction horizon, control horizon, solver iteration limits, and constraint penalty tuning; these parameters often require re-adjustment as the sampling time or operating regime changes. Most crucially however, poor initialisation or an infeasible problem can cause solver failure in exactly the high-demand conditions where reliable control matters most; a critical issue for open-loop-unstable systems like this hexacopter.

SDRE-family methods occupy the middle ground: they keep the LQ structure (and its tooling, intuition, and predictable compute) while re-'linearising' *exactly* — the SDC factorisation $f(x) = A(x)x$ is an algebraic identity, not a Taylor truncation, at every step. The contribution here is a tracking formulation of that idea that is receding-horizon, preview-aware, numerically robust over long runs, and deliberately engineered around a fixed sequence of dense linear-algebra operations, yielding *consistent per-step timing* and *graceful behaviour at the edges* (end of reference, preview window too short, loss of preview, solver early termination).

The proposed SDDRE-based approach avoids NLMPC's optimisation-tuning loop by replacing online nonlinear programming with a fixed sequence of Riccati-based computations. With a fixed Newton–Kleinman iteration budget and bounded preview recursion, the computational workload is known a priori, providing consistent execution timing while retaining near-optimal nonlinear tracking capability, and enabling fast sampling rates.

---

## The controller

### Problem setup

State $x = [\,p_I;\ v_b;\ \Theta;\ \omega\,] \in \mathbb{R}^{12}$ (inertial position, body-frame velocity, attitude, body rates); input $u = \delta\omega \in \mathbb{R}^6$, rotor-speed deviations from hover. Full nonlinear rigid-body dynamics with a quadratic rotor thrust/torque map; the simulation truth model integrates these with `ode45` under ZOH inputs and actuator saturation.

Given a reference $\{r_k\}_{k=0}^{N}$ for outputs $y = Cx$ (position, yaw, and roll/pitch rates in the benchmark), minimise

$$
J \;=\; \sum_{k=0}^{N-1} \Big[ (Cx_k - r_k)^T Q_y (Cx_k - r_k) + u_k^T R\, u_k \Big]
\;+\; (Cx_N - r_N)^T Q_{yf} (Cx_N - r_N).
$$

Weights are Bryson-scaled from physically meaningful allowable deviations and the vehicle's RPM headroom.

### SDC factorisation

State Dependent Coefficient (SDC) factorising: the nonlinear dynamics are rewritten **exactly** as $\dot x = A(x)\,x + B(u)\,u$:

- $B(u)$ is exact, not a Jacobian: the thrust increment obeys
$k_F\big((\bar\omega + \delta)^2 - \bar\omega^2\big) = k_F(2\bar\omega + \delta)\,\delta$, so the input matrix carries
$(2\bar\omega_i + \delta u_i)$ evaluated at the previous input.
- The gravity residual $g\,(e_3 - C_{bI}e_3)$ is factorised smoothly through zero attitude using 
$\mathrm{sinc}$ and versine forms, e.g.
$A_{4,8} = g\,\mathrm{sinc}(\theta)$,
$A_{6,7} = g\,\frac{1-\cos\phi}{\phi}$ — no singularities, exact at every operating point.
- Attitude kinematics and gyroscopic coupling are factorised similarly. Four attitude parameterisations are implemented (Euler, quaternion, MRP, and exponential coordinates / FRA); the benchmark uses the **FRA** (*Finite Rotation Angle*) factorisation, which degrades most gracefully at the performance limit and requires no attitude-weight rescaling relative to the Euler design.

SDC factorisations are non-unique; the choices above were selected empirically and each is documented at its definition site.

### Per-step algorithm

At step $k$, with measured state $x_k$ and previous input $u_{k-1}$:

1. **Factorise & discretise.** $A_k = A(x_k)$, $B_k = B(u_{k-1})$; ZOH-discretise at $T_s$.
2. **Riccati solve (warm-started).** Solve the DARE
$P = A^TPA - A^TPB(R + B^TPB)^{-1}B^TPA + C^TQ_yC$ for $P_{ss}$ using Cascade Newton–Kleinman (C-NK): initialise from the previous step's $P_{ss}$, then iterate
$K \leftarrow (R + B^TPB)^{-1}B^TPA$, $\;\;P \leftarrow \mathrm{dlyap}\big((A - BK)^T,\ Q + K^TRK\big)$
(quadratically convergent, see [3-5]).
Because consecutive SDC models differ little, a few iterations typically suffice.
Set $K_{ss}$, $A_{cl} = A - BK_{ss}$.
3. **Preview feedforward.** Over a receding window of $M$ steps, run the costate recursion
$v \leftarrow A_{cl}^T v + C^TQ_y\,r_{k+j}$ backwards from the seed
$$v_{k+M} \;=\; (I - A_{cl}^T)^{-1} C^T Q_y\, r_{k+M}$$ (see *Graceful degradation*). The control is then: $$u_k = -K_{ss}\,x_k \;+\; (R + B^TP_{ss}B)^{-1}B^T v_{k+1}.$$
4. **Terminal regime.** Once the window reaches the end of the trajectory ($k + M \ge N$), the controller switches to the **full finite-horizon difference-Riccati recursion** from $C^TQ_{yf}C$ — exactly the frozen-SDC LQT solution — so the gain ramps correctly into the terminal cost. (Disable via `UseFullFiniteHorizonMPCAtTerminal=false` if strictly constant solve time is preferred over exact terminal handling.)

> **Terminology.** "SDDRE" follows Khamis & Naidu's finite-horizon SDRE tracking framework [2], transposed to discrete time. In the standing regime the controller solves a state-dependent *algebraic* DARE with a preview costate; the difference-Riccati recursion proper is engaged in the terminal regime (and available throughout via `AlwaysUseFullFiniteHorizonMPC`).

---

## Results

### Feasible maneuver envelope

The main test trajectory is a 3D maneuver along a closed path of ≈30 m arc length[^arc_length], with a figure-eight lateral component:
$$r(t) = 4\,[\sin\omega t,\ -\cos\omega t\,\sin\omega t,\ \cos\omega t], \quad \omega = 2\pi/t_f$$

swept over duration $t_f \in [4, 20]$ s. Identical cost matrices for both controllers.

[^arc_length]: Path length: 16√2·E(1/√2) ≈ 30.562 m analytically, while empirically a well-tracked trajectory measured ≈30.564 m.

<p align="center">
  <img src="docs/figures/errors_and_efforts.png" width="600"
       alt="Position RMSE and control effort vs maneuver time, SDDRE vs Linear"/>
</p>

The two controllers are near-indistinguishable while the small-angle assumption holds; the gap opens exactly where the linear model stops being true.

On amplitude sweeps of a circle-like reference at 50 Hz, the linear design fails near amplitude 1.7m while SDDRE fails near
3.1m (Euler SDC), 
3.4m (quaternion SDC), 
3.6m (MRP SDC) and 
3.7m (FRA SDC).

### Consistent timing

Per-step compute is a fixed pipeline: an SDC evaluation, discretisation, warm-started DARE solve, and a fixed-length ($M$-step) costate recursion. Total simulation time for the benchmark maneuver, by DARE strategy:

| DARE solver | $T_s = 1/200$ (desktop PC) |
|---|---|
| `dlqr` (cold)  | 7.25 s |
| `idare` (cold) | 5.80 s |
| **Cascade Newton–Kleinman (warm)** | **4.15 s** |
| (no DARE solves — floor) | 3.70 s |

Warm-starting removes ~80–85 % of the Riccati cost; what remains is dominated by the (fixed-cost) preview recursion and other overhead.
> Note that these figures record full simulation times, with `ode45` overhead being the bulk of this, hence the no-DARE-solves floor was included as an important reference. A more precise benchmark isolating the controller compute time and comparing across LQT/SDDRE/NLMPC schemes, including timing variability, will be placed here in future.


Two properties make this more than an average-case speedup:

- **Bounded iterations.** Every Newton–Kleinman iterate from a stabilising seed is itself a   stabilising solution estimate (Hewer [5]), and the iterates decrease monotonically to $P_{ss}$. An iteration cap therefore yields a *hard per-step compute bound* while the resulting gain remains stabilising for the frozen model — early termination is suboptimal but safe.
- **No structural variability.** There is no active set, no line search, no initialisation-dependent NLP path. The operation count has a hard, state-independent upper bound (capped NK iterations, fixed-length preview recursion) which can be further bounded or fixed so it's completely known beforehand, ensuring timing consistency and robustness. No difficulty starting cold, convergence issues, or tuning/validation of hyperparameters. There is no NLP to solve — the per-step workload is replaced by deterministic linear algebra, easily deployed and verified on a microcontroller

### Graceful degradation

Receding-horizon controllers must behave well when preview information is 'limited' (when the reference is only known for a short window ahead), or when it runs out (at the terminal condition). Here that property is obtained by construction: the preview recursion is seeded with the fixed point of its own constant-tail dynamics, $(I - A_{cl}^T)^{-1}C^TQ_y\,r_{k+M}$, so that as the window slides off the end of the reference (or as $M \to \infty$, or $M \to$ small), the feedforward collapses continuously onto the infinite-horizon preview solution. This is a structural property, not dependent on solver convergence or feasibility.

In NLMPC, a shrinking or insufficient horizon can render the NLP infeasible, producing no control action precisely at the trajectory boundary; here, the controller always returns a stabilising gain and a feedforward that smoothly interpolates between the finite-horizon and infinite-horizon solutions. The naive seed ($C^TQ_y\,r_{k+M}$​) lacks this fixed-point consistency: it implicitly assumes the costate is zero beyond the window, which is only true at a genuine trajectory endpoint. When the window slides past the end of the reference, or when the preview horizon is short relative to the closed-loop settling time, this mismatch injects a transient into the feedforward at every step.


## Usage

### Requirements

- MATLAB R2021a or newer (name–value function-call syntax)
- Control System Toolbox (`ss`, `c2d`, `idare`, `dlyap`)
- ... no Optimization Toolbox, YALMIP, or external solvers required[^mpc_toolbox_req]

[^mpc_toolbox_req]: Note that future versions may require MATLAB's MPC Toolbox in order to run NLMPC tracking comparisons (SDDRE vs NLMPC).

### Quickstart

See `demo_sddre.m` for the main work. Here you can run LQT and the SDDRE controller over the sample trajectory, inspecting the results — feel free to play around with sample rate `Ts`, cost matrices `Q` and `R`, `PreviewHorizon`, and the actual reference trajectory to be tracked. MATLAB's NLMPC control will also be exposed via simple wrapper for a full comparison in a future release. Example usage:

```matlab
% Plant, rates, cost (see demo_sddre.m for the full setup)
load_copter_params;   qp.Ts = 1/50;
qp.nominal_omegas = compute_omega_bar(qp, ones(1,qp.n_rotors), 'least_squares');
% ... build Q, R, C, Qy from allowable deviations (demo shows the Bryson scaling) ...

% Reference: figure-eight-style 3D maneuver, completed over t_f seconds, sampled on the control grid
tspan = [0 8.5];
w = 2*pi/tspan(2);
x_ref_fcn = @(t) 4*[sin(t*w), -cos(t*w)*sin(t*w), cos(t*w), 0 0 0, 0 0 0, 0 0 0]';
tgrid = 0:qp.Ts:tspan(2);
r_ = C * cell2mat(arrayfun(x_ref_fcn, tgrid, 'UniformOutput', false));

% SDDRE controller (FRA attitude SDC, 2s preview)
ctrl_fcn = @(t,x,k,u_prev) compute_u_SDDRE_v3(t, x, k, u_prev, r_, C, Qy, R, 100*Qy, ...
    tspan(2), qp, PreviewHorizon=2.0, SDC_A_function=@get_A_matrix_SDRE_FRAAttitude);

% Simulate against the full nonlinear model (ode45 truth, ZOH inputs, actuator limits)
[t, X, U] = run_sim(qp, tspan, x0, @(t) ones(1,6), ctrl_fcn, AttitudeRepresentation='fra');
do_plots
```

### Repository layout

```
demo_sddre.m                 Minimal end-to-end demo (the quickstart, runnable)
benchmark_sweep.m            Reproduces the maneuver-time sweep figure
compute_u_SDDRE_v3.m         The controller
iterative_dare.m             Warm-started DARE solver (Riccati iteration / Newton–Kleinman)
get_A_matrix_SDRE_*.m        SDC factorisations (Euler / Quaternion / MRP / FRA)
get_B_matrix_SDRE.m          Exact input-matrix SDC
solve_LQT.m                  Finite-horizon LQT (the linear baseline + terminal-regime math)
run_sim.m, nonlinear_dynamics.m, load_copter_params.m, compute_omega_bar.m
                             Simulation truth: full 6DOF model, 4 attitude parameterisations
docs/figures/                Figures used in this README + benchmark data (.mat)
```

To reproduce the headline figure: run `benchmark_sweep.m` (sweeps $t_f$, runs both controllers, saves `results_*.mat`), then the plotting cell at its end. Fair warning: the full sweep is ~90 `ode45` simulations and takes a while to compute.


## How this solution was arrived at

This controller is the survivor of a broader exploration of the tracking problem, roughly in order:

- **Linear LQT** with an implicit-reference/feedforward split — excellent inside the small-angle regime, and the structural backbone everything else reuses.
- **Robust linear tracking** against the uncertainty polytope: brute-force minimax Riccati recursions, per-timestep LMI synthesis ($H_2$, $H_\infty$, guaranteed-cost; quadratic-stability vs parameter-dependent Lyapunov), trajectory-driven uncertainty-set construction, μ-analysis and disk margins. Effective for genuine uncertainty — but robustness to *known* effects from nonlinearity from large attitudes is the wrong tool: the linear model isn't uncertain there, it's wrong in a known way. While the practical payoff of a robust method working is huge (a fixed feedback gain Kr across the whole trajectory, and a simple LQT-like recursion for the feedforward), the conservatism it introduces proves too problematic.
- **Iterative methods**: iLQR - not suitable for online control. Marginally improves tracking, but computationally heavy and trades away the properties that matter for flight code (convergence reliability, tuning burden, and efficient, predictable compute).
- **Nonlinear MPC** (`nlmpc`, and `nlmpcMultistage`): the accuracy ceiling, and a valuable benchmark — but with the online-NLP timing characteristics, compute cost, hyperparameter sensitivity, and convergence issues discussed above.
- **SDDRE**, which keeps the LQ machinery, considers the exact nonlinear dynamics, previews the reference, and budgets its compute whilst maintaining reliability — the sort of approach this project was looking for.


## Limitations & honest caveats

- All results are simulation-based with the setup previously described — aerodynamic disturbances, energy dissipation, and other more sophisticated dynamics are neglected.
- No actuator *rate* constraints are assumed or modelled. System parameters (moments of inertia, thrust coefficients) are assumed to be known without error.
- Full state measurement is assumed available and free from error — assumes no estimator dynamics or noise in the state vector.
- SDRE-family methods carry no *global* optimality or stability guarantee — pointwise stabilisability of the SDC pair is assumed (and monitored), and the SDC factorisation choice is a genuine design degree of freedom. Different parametrisations are mathematically valid and pointwise-equivalent, yet will yield slightly different behaviour — requires empirical testing and a bit of finesse (Cloutier referred to this aspect as an 'art' [6]).
- The DARE warm start assumes model continuity between steps; a discontinuous state jump falls back to a cold `idare` solve (handled, but that step is slower).

## Next steps

- **acados NLMPC comparison** — a proper real-time-iteration SQP baseline, to put honest numbers on the compute-performance-reliability frontier this controller claims to sit on, as well as stress-testing NLMPC's practical performance (tendency to get stuck in local minima, difficulty converging or finding initial solution, infeasibility – empirically observed in MATLAB's NLMPC implementations, but these are not industry-standard NLMPC solvers – would need more solid grounding for that claim).
- Hardware-oriented port: the per-step pipeline (SDC → ZOH → capped NK → fixed recursion) is embedded-friendly by design; a C implementation with static allocation is the natural test of the timing claims.

## References

[1] B. D. O. Anderson, J. B. Moore, *Optimal Control: Linear Quadratic Methods*, Prentice Hall, 1990.

[2] A. Khamis, D. S. Naidu, *Nonlinear optimal tracking using finite-horizon State Dependent Riccati Equation (SDRE)*, 2014. doi:10.1016/j.isatra.2014.06.006

[3] L. Saluzzi, M. Strazzullo, *Dynamical Low-Rank Approximation Strategies for Nonlinear Feedback Control Problems*, 2025. arXiv:2501.07439.

[4] L. Saluzzi, *The State-Dependent Riccati Equation in Nonlinear Optimal Control: Analysis, Error Estimation and Numerical Approximation*, 2025. arXiv:2503.01587.

[5] G. Hewer, *An iterative technique for the computation of the steady state gains for the discrete optimal regulator*, 1971 (Discrete Newton–Kleinman, the properties underwriting the iteration-cap argument.). doi:10.1109/TAC.1971.1099755

[6] J. R. Cloutier, D. T. Stansbery, *The capabilities and art of state-dependent Riccati equation-based design*, 2002. doi:10.1109/ACC.2002.1024785

[7] T. Çimen, *State-Dependent Riccati Equation (SDRE) Control: A Survey*, 2008. doi:10.3182/20080706-5-KR-1001.00635

[8] S. W. Hur, S. H. Lee, C.‑J. Kim, *Effects of Attitude Parameterization Methods on Attitude Controller Performance*, 2020. doi:10.1007/s42405-020-00286-3