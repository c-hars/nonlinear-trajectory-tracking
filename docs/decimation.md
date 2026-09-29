## Motivation

The preview control law

$$
u_k = -K x_k + K_v v_{k+1}
$$

admits a clean split that lets it scale well to high loop rates.

Consider kHz-rate control of the hexacopter. Computing and applying the feedback $-K x$ at 1000Hz has a clear advantage over say 20Hz control, but the reference preview $K_v v_{k+1}$ does not benefit in the same way; it benefits from the finer loop rate only insofar as the continuous-time dynamics are better captured by the 1000Hz model over the 20Hz model. For the hexacopter, this is a diminishing return. There's little difference between the quality of a 20Hz reference preview and a 1000Hz reference preview, since its dynamics and modes are already low-frequency enough. The net effect is largely just a 50x increase in compute, rather than improved tracking, and this prevents high loop rate control.

"Decimation" of the reference preview is the solution to this. At each control step, update and apply the gain $K$ for feedback, but compute the reference preview $K_v v_{k+1}$ over a *decimated grid*; group multiple fine steps into coarser "blocks", and compute the recursion over these. Effectively: feedback that scales with loop rate, preview that doesn't – the benefits of high loop rate feedback control without the largely-redundant overhead.

A similar decimation scheme applies to the shrinking horizon / MPC mode, where the full Riccati recursion is needed (the time-varying cost-to-go/gain, as well as the costate).

## Glossary: "decimated recursions" and multi-rate MPC

Terminology, useful for understanding the functions in the repo and the math that follows.

- **Fine step:** a prediction step of duration `Ts`.
- **Coarse step:** a prediction step covering more time than a fine step (duration > `Ts`); lower model fidelity, faster compute. A coarse step has duration `d*Ts`, where `d` is the decimation factor; `d` consecutive fine steps are grouped and approximated as one coarse step.
- **Block:** a discrete-time step over `d*Ts` seconds; 'block' is shorthand for this chunk of time.
- **Coarse matrices:** `(A_d, B_d, Qy_d, R_d)` – the system matrices `(A, B, Q, R)` after applying a decimation factor. They represent the same† continuous-time dynamics and OCP formulation `(Ac, Bc, Qy_c, R_c)`, but discretised at a different (lower-resolution) sampling rate. Also referred to as `d`-step matrices.
- **Fine recursion:** a backward recursion processing one fine step at a time.
- **Coarse recursion:** a backward recursion processing one coarse step at a time.
- **Tier:** a contiguous sequence of steps sharing the same decimation factor.
  - **Far-horizon tier:** the final tier – farthest from the current time, and uses the largest decimation factor. Absorbs whatever horizon remains after the near-horizon tiers.
  - **Near-horizon tier:** any tier other than the far-horizon tier. Each holds `StepsPerTier` steps; each step spans a block of `d` fine steps.
- **Remainder fine steps:** fine steps that don't fit cleanly into the tier schedule; the prediction horizon of `M` steps might have a remainder. These are not discarded, instead they're placed nearest to the current time and processed by the fine recursion.
- **Preview grid:** the full arrangement across the preview horizon. Forward in time: remainder fine steps → near-horizon tiers → far-horizon tier → preview endpoint. The backward recursion processes this in reverse.

> (†) not exactly the same: some approximation also occurs, described further on.

# Summary of formulas

**Notation**

$(A,B)$ are the system's discrete-time matrices at the nominal $T_s$ – exact discretisations of $(A_c,B_c)$, the continuous-time matrices.

$K$ refers to the gain associated with $(A,B,Q,R)$

A subscript of $d$ refers to the decimated matrices – $d$-step system matrices.

Shorthand, repeated throughout:

$$
Q = C^T Q_y C, \qquad q_j = C^T Q_y r_{k+j}
$$

$$
\bar{q} = C^T Q_y \bar{r}, \qquad \bar{r} = \frac{1}{d} \sum_{i=0}^{d-1} r_{k+j+i}
$$


Note that the costate/cost-to-go seeds are not stated here for brevity.


## Preview control

Set $P$, $K$, $A_{cl}$ via the DARE.

Both preview control cases use the control input:

$$
u_k = -K x_k + (R + B^T P B)^{-1} B^T v
$$

where $v$ is the costate, computed via backward recursion over the preview horizon – in fine steps, coarse steps, or some mixture of both; however the grid structure is specified.

### Standard ("fine") recursion

The recursion is:
$$
v \leftarrow F v + q_j
$$

where $F = A_{cl}^T$.

$\square$

### $d$-step recursion ("coarse" recursion)

Assumption: the fine-rate feedback $u = -K x$ *acts at every step* throughout the block, i.e.
$$
u_{k+j+i} = -K x_{k+j+i}
$$
for $i = 0, \dots, d-1$ and some $j \geq 0$.

The exact $d$-step formula is then:
$$
v \leftarrow F^d v + \sum_{i=0}^{d-1} F^i q_{j+i}
$$
With $r$ held at its block average, this becomes:
$$
v \leftarrow F_d v + G_d \bar{q}
$$

where
$F_d = F^d$ and $ G_d = \sum_{i=0}^{d-1} F^i $.

$\square$

## MPC

### Standard ("fine") recursion

$$
\begin{aligned}
K &= (R + B^T P B)^{-1} B^T P A \\
A_{cl} &= A - BK \\
P &\leftarrow Q + K^T R K + A_{cl}^T P A_{cl} \\
v &\leftarrow A_{cl}^T v + q_j \\
\end{aligned}
$$
$$
u_k = -(R + B^T P B)^{-1} (B^T P A x_k - B^T v)
$$

$\square$


### $d$-step recursion ("coarse" recursion)

Assumption: coarse-rate feedback $u = -K_d x$ is *held* across the block, i.e.
$$
u_{k+j+i} = u_{k+j} = -K_d x_{k+j}
$$
for $i = 0, \dots, d-1$ and some $j \geq 0$. The exact $d$-step formula for recursion is then:

$$
\begin{aligned}
K_d &= (R_d + B_d^T P B_d)^{-1} (B_d^T P A_d + N_d^T) \\
A_{cl,d} &= A_d - B_d K_d \\
P &\leftarrow Q_d + K_d^T R_d K_d - N_d K_d - K_d^T N_d^T + A_{cl,d}^T P A_{cl,d} \\
v &\leftarrow A_{cl,d}^T v + \sum_{i=0}^{d-1} (A_i - B_i K_d)^T q_{j+i}
\end{aligned}
$$

where
$$
\begin{aligned}
A_i &= A^i, & A_d &= A_i|_{i=d} \\
B_i &= (I + A + \dots + A^{i-1}) B, & B_d &= B_i|_{i=d} \\
\end{aligned}
$$
$$
Q_d = \sum_{i=0}^{d-1} A_i^T Q A_i, \qquad N_d = \sum_{i=0}^{d-1} A_i^T Q B_i, \qquad
R_d = d R + \sum_{i=0}^{d-1} B_i^T Q B_i
$$

With $r$ held, the costate update simplifies to:
$$
v \leftarrow A_{cl,d}^T v + (G_d - H_d K_d)^T \bar{q}
$$

where
$G_d = \sum_{i=0}^{d-1} A_i$ and $H_d = \sum_{i=0}^{d-1} B_i$. Finally, with the held-$x$ assumption (stage cost applies as if there's no fine-rate dynamics within a block, i.e. only applies at the block rate), the recursion further simplifies:
$Q_d = dQ, \ R_d = dR, \ G_d = dI, \ N_d = H_d = 0$, and with these it becomes the familiar form:
$$
\begin{aligned}
K_d &= (R_d + B_d^T P B_d)^{-1} B_d^T P A_d \\
A_{cl,d} &= A_d - B_d K_d \\
P &\leftarrow Q_d + K_d^T R_d K_d + A_{cl,d}^T P A_{cl,d} \\
v &\leftarrow A_{cl,d}^T v + d \, \bar{q} \\
\end{aligned}
$$

Note that if the first block is coarse ($j = 0$), the control input would then be:
$$
u_k = -(R_d + B_d^T P B_d)^{-1} (B_d^T P A_d x_k - B_d^T v)
$$

$\square$





# Supporting math.

The system (undecimated) is:
$$
x_{k+1} = A x_k + B u_k
$$
$$
J = \sum_{k=0}^N (Cx_k-r_k)^T Q_y (Cx_k-r_k) + u_k^T R u_k 
$$

---

### $d$-step preview recursion

After seeding, the standard recursion update is
$$ v \leftarrow A_{cl}^T v + C^T Q_y r_{k+j} $$
where $A_{cl} = A-BK$ and $K=K_{ss}$ via the DARE.

More concisely, writing $F = A_{cl}^T$ and $q_j = C^T Q_y r_{k+j}$ for brevity, we'll write this as:
$$ v \leftarrow F v + q_j $$

Now, consider a d-step recursion

$$ v \leftarrow F v + q_{j+d-1} $$
$$ \vdots $$
$$ v \leftarrow F v + q_{j} $$

i.e.

$$
\begin{aligned}
\text{Step 1} \rightarrow v &= Fv + q_{j+d-1} \\
\text{Step 2} \rightarrow v &= F(Fv + q_{j+d-1}) + q_{j+d-2} \\
  &= F^2v + Fq_{j+d-1} + q_{j+d-2} \\
\text{Step 3} \rightarrow v &= F(F^2v + Fq_{j+d-1} + q_{j+d-2}) + q_{j+d-3} \\
&= F^3v + F^2q_{j+d-1} + Fq_{j+d-2} + q_{j+d-3} \\
&\vdots \\
\text{Step }d \rightarrow v &= F^dv + F^{d-1}q_{j+d-1} + ... + Fq_{j+1} + q_j \\
\end{aligned}
$$

So the sequence over a $d$-step block collapses to, equivalently, one update involving $d$ matvecs:
$$
v \leftarrow F^d v + \sum_{i=0}^{d-1} F^i q_{j+i}
$$

Finally, let $r$ be held fixed over the block, $r=\bar{r}$. Then the $d$ fine steps collapse into one update involving $2$ matvecs:
$$
v \leftarrow F_d v + G_d \bar{q}
$$
where $F_d = F^d$, $ G_d = \sum_{i=0}^{d-1} F^i $, and $\bar{q} = C^T Q_y \bar{r}$.

---

### $d$-step matrices $(A_d,B_d)$

Consider the dynamics' evolution over $d$ steps:

$$
\begin{aligned}
x_{k+1} &= A x_k + B u_k \\
x_{k+2} &= A x_{k+1} + B u_{k+1} \\
&= A (A x_k + B u_k) + B u_{k+1} \\
&= A^2 x_k + AB u_k + B u_{k+1} \\
x_{k+3} &= A x_{k+2} + B u_{k+2} \\
&= A^3 x_k + A^2B u_k + AB u_{k+1} + B u_{k+2} \\
x_{k+4} &= A^4 x_k + A^3B u_k + A^2B u_{k+1} + AB u_{k+2} + B u_{k+3} \\
&\vdots \\
x_{k+d} &= A^d x_k + \sum_{i=0}^{d-1} A^i B u_{k+d-1-i}
\end{aligned}
$$
Assume the input is held (ZOH). Then $u_k = u_{k+1} = ...$, so:
$$ x_{k+d} = \underbrace{A^d}_{:= A_d} x_k + \underbrace{(I+A+A^2+...+A^{d-1}) B}_{:= B_d} u_k $$

---

### $d$-step Riccati recursion

Recall that the original OCP is defined with respect to
$$
x_{k+1} = A x_k + B u_k
$$
$$
J = \sum_{k=0}^N (Cx_k-r_k)^T Q_y (Cx_k-r_k) + u_k^T R u_k 
$$

and that under the held-input assumption (*move blocking*),
$$
x_{k+i} = A_i x_k + B_i u_k
$$
where $ A_i = A^i $ and $ B_i = (I+A+A^2+...+A^{i-1}) B $.

Let's consider the accumulated cost over $d$ steps.

At step $k$ we have
$$
J_k = x_k^T C^T Q_y C x_k - 2 x_k^T C^T Q_y r_k + r_k^T Q_y r_k + u_k^T R u_k
$$

This breaks down into, respectively: *state* cost, *costate* cost, a constant (optimal input is invariant w.r.t this), and *input* cost. We split the analysis into the first two cases.

---

The accumulated ***state*** cost over $d$ stages is
$$
\sum_{i=0}^{d-1} x_{k+i}^T Q x_{k+i} \qquad (†)
$$
where $Q$ = $C^T Q_y C$. And under the held-input assumption, this becomes
$$
\sum_{i=0}^{d-1} (A_i x_k + B_i u_k)^T Q (A_i x_k + B_i u_k)
$$
$$
= x_k^T Q_d x_k + 2 x_k^T N_d u_k + u_k^T S_d u_k
$$
where
$$
\begin{aligned}
Q_d &= \sum_{i=0}^{d-1} A_i^T Q A_i \\
N_d &= \sum_{i=0}^{d-1} A_i^T Q B_i  \\
S_d &= \sum_{i=0}^{d-1} B_i^T Q B_i \\
\end{aligned}
$$
This is the *exact* state cost accumulated over $d$ steps.

If we then introduce a held-x assumption (along with the held-u assumption), then we get this simplified; there are no dynamics within a block, i.e. $A_i \rightarrow I$ and $B_i \rightarrow 0 $, so we get only $d$ repetitions of $Q$ and (†) becomes just
$$ x_k^T Q_d x_k $$
where $ Q_d = dQ $. The cross term $N_d$ and the input term $S_d$ both vanish.

The held-x assumption is used in the implementation for a few reasons. Apart from maintaining simplicity, the held-x assumption also maintains the spirit of decimation: far in the prediction horizon, trust the prediction model less. That includes the state evolution model (held-x assumption), not just move blocking (held-u assumption). It's a cleaner design that also tends to yield slightly better tracking (and slightly faster compute).

---

The accumulated ***costate*** cost over $d$ stages is
$$
-2 \sum_{i=0}^{d-1} x_{k+i}^T C^T Q_y r_{k+i} \qquad (‡)
$$

Similarly, under the held-u assumption the sum becomes:
$$
\sum_{i=0}^{d-1} (A_i x_k + B_i u_k)^T C^T Q_y r_{k+i}
$$

Let $q_{k+i} = C^T Q_y r_{k+i}$ for brevity. Note that $u_k$ is fixed here under the held-u assumption: $u_k = - K_d x_k$, where $K_d$ is the coarse gain, so this becomes:
$$
= x_k^T \sum_{i=0}^{d-1} (A_i - B_i K_d)^T q_{k+i}
$$

which is the *exact* costate cost accumulated over $d$ steps.

The costate from the next block, $v$, enters the cost as $-2 v^T x_{k+d}$, and with $x_{k+d} = A_d x_k + B_d u_k = A_{cl,d} x_k$ it contributes $A_{cl,d}^T v$. The $d$-step costate update is therefore:
$$
v \leftarrow A_{cl,d}^T v + \sum_{i=0}^{d-1} (A_i - B_i K_d)^T q_{k+i}
$$

If we also hold $r$ at its block average, $r_{k+i} = \bar{r}$, so $q_{k+i} = \bar{q} = C^T Q_y \bar{r}$, then the sum factors cleanly.

So with $r$ held, the $d$-step costate update becomes:
$$
v \leftarrow A_{cl,d}^T v + (G_d - H_d K_d)^T \bar{q}
$$
where $G_d = \sum_{i=0}^{d-1} A_i$ and $H_d = \sum_{i=0}^{d-1} B_i$.

This is the *exact* expression for the optimal $d$-step costate update, accounting for steps $k,k+1,...,k+d-1$, under held $u_k=-K_d x_k$ and held $r=\bar{r}$.

--

Now, under an additional held-x assumption, $A_i \rightarrow I$ and $B_i \rightarrow 0$ within the block, and so the sum becomes:
$$
x_k^T \sum_{i=0}^{d-1} C^T Q_y r_{k+i} = x_k^T C^T (d Q_y) \bar{r}
$$
so the $d$-step costate update is then simply:
$$
v \leftarrow A_{cl,d}^T v + C^T (d Q_y) \bar{r}
$$
Note that under held-x, the block averaging of $r$ becomes *exact* rather than an assumption: with the state held, only the sum of the references matters. This is another reason why the held-x assumption was chosen; it opens up a clean way of decimating $r$ without full information loss; block-averaging instead of downsampling.

## Appendix: $d$-step Riccati recursion: full derivation (including within-block feedforward)

The above math writes the input as $u = -K_d x$, neglecting the feedforward terms to keep things clean – they don't affect the $P$ and $v$ recursions, since their cross terms cancel by optimality of $K_d$. The below is a more rigorous first-principles derivation that includes them.

*Attribution: Co-authored by Claude Opus 5.5. Edited and reviewed by a human's large language model.*

**Setup.** Consider one block starting at offset $j$, with the input held across it: $u_{k+j+i} = u$ for $i = 0, \dots, d-1$. Write $x = x_{k+j}$ for the state at the block start. The fine stage cost is
$$
\begin{aligned}
\ell_{k+j+i} &= (C x_{k+j+i} - r_{k+j+i})^T Q_y (C x_{k+j+i} - r_{k+j+i}) + u^T R u \\
&= x_{k+j+i}^T Q x_{k+j+i} - 2 q_{j+i}^T x_{k+j+i} + u^T R u + \text{const}
\end{aligned}
$$
and the cost-to-go at the end of the block (offset $j+d$) is, up to a constant,
$$
V^+(x^+) = x^{+T} P x^+ - 2 v^T x^+
$$
Terms constant in $(x, u)$ don't affect the gain or costate, so they're dropped throughout.

**Within-block states.** Applying the fine dynamics $i$ times with $u$ held:
$$
x_{k+j+i} = A_i x + B_i u, \qquad x^+ = A_d x + B_d u
$$
with $A_i = A^i$, $B_i = (I + A + \dots + A^{i-1}) B$, and $A_0 = I$, $B_0 = 0$.

**Block cost.** Take a single fine stage $i$ and substitute $x_{k+j+i} = A_i x + B_i u$ into each of its three terms.

The quadratic state term expands as
$$
(A_i x + B_i u)^T Q (A_i x + B_i u) = x^T (A_i^T Q A_i) x + 2 x^T (A_i^T Q B_i) u + u^T (B_i^T Q B_i) u
$$
The linear (reference) term becomes
$$
-2 q_{j+i}^T (A_i x + B_i u) = -2 (A_i^T q_{j+i})^T x - 2 (B_i^T q_{j+i})^T u
$$
and the input term $u^T R u$ is unchanged, since $u$ is held.

Summing over $i = 0, \dots, d-1$ and collecting like terms:
$$
\sum_{i=0}^{d-1} \ell_{k+j+i} = x^T \underbrace{\left( \sum_i A_i^T Q A_i \right)}_{Q_d} x + 2 x^T \underbrace{\left( \sum_i A_i^T Q B_i \right)}_{N_d} u + u^T \underbrace{\left( dR + \sum_i B_i^T Q B_i \right)}_{R_d} u - 2 \underbrace{\left( \sum_i A_i^T q_{j+i} \right)^T}_{s_x^T} x - 2 \underbrace{\left( \sum_i B_i^T q_{j+i} \right)^T}_{s_u^T} u
$$
So $Q_d$, $N_d$, $R_d$ collect the quadratic terms (as found previously), and the reference terms collect into
$$
s_x = \sum_{i=0}^{d-1} A_i^T q_{j+i}, \qquad s_u = \sum_{i=0}^{d-1} B_i^T q_{j+i}
$$

**Total cost from the block start.** The total cost from the block start is the block cost plus the cost-to-go from the end of the block:
$$
J(x, u) = \sum_{i=0}^{d-1} \ell_{k+j+i} + V^+(x^+), \qquad x^+ = A_d x + B_d u
$$
Expanding $V^+$ we have:
$$
V^+(A_d x + B_d u) = x^T (A_d^T P A_d) x + 2 x^T (A_d^T P B_d) u + u^T (B_d^T P B_d) u - 2 (A_d^T v)^T x - 2 (B_d^T v)^T u
$$
Adding this to the block cost and collecting like terms:
$$
J(x, u) = x^T (Q_d + A_d^T P A_d) x + 2 x^T (N_d + A_d^T P B_d) u + u^T (R_d + B_d^T P B_d) u - 2 (s_x + A_d^T v)^T x - 2 (s_u + B_d^T v)^T u
$$

**Optimal input.** Differentiating $J$ with respect to $u$:
$$
\frac{\partial J}{\partial u} = 2 (N_d + A_d^T P B_d)^T x + 2 (R_d + B_d^T P B_d) u - 2 (s_u + B_d^T v) = 0
$$
Rearranging:
$$
(R_d + B_d^T P B_d) u = -(B_d^T P A_d + N_d^T) x + (B_d^T v + s_u)
$$
i.e.
$$
u^* = -K_d x + (R_d + B_d^T P B_d)^{-1} (B_d^T v + s_u)
$$
where
$
K_d = (R_d + B_d^T P B_d)^{-1} (B_d^T P A_d + N_d^T)
$.

The feedforward has two parts. $B_d^T v$ accounts for the reference after the block: $B_d$ is how $u$ moves the state by the block's end, and $v$ carries the reference information from there on. $s_u$ is the *within-block feedforward*: it accounts for the reference during the block, since $B_i$ is how the held $u$ moves the state by each intermediate stage $i$.

**Cost-to-go update.** The cost-to-go at the block start is $\min_u J(x, u)$, written as a function of $x$ alone. Start from $J$ as derived above:
$$
J(x, u) = x^T (Q_d + A_d^T P A_d) x + 2 x^T (N_d + A_d^T P B_d) u + u^T (R_d + B_d^T P B_d) u - 2 (s_x + A_d^T v)^T x - 2 (s_u + B_d^T v)^T u
$$
For brevity, write $L = N_d + A_d^T P B_d$ and $w = s_u + B_d^T v$. The two $u$-linear terms then combine, $2 x^T L u - 2 w^T u = 2 (L^T x - w)^T u$, so $J$ splits into:
$$
J(x, u) = \underbrace{x^T (Q_d + A_d^T P A_d) x - 2 (s_x + A_d^T v)^T x}_{\text{no } u} + \underbrace{u^T (R_d + B_d^T P B_d) u + 2 (L^T x - w)^T u}_{u\text{-dependent}}
$$

Substituting $u^* = -(R_d + B_d^T P B_d)^{-1} (L^T x - w)$ into the $u$-dependent part of $J$, then using $(R_d + B_d^T P B_d)^{-1} L^T = K_d$:
$$
\begin{aligned}
&u^{*T} (R_d + B_d^T P B_d) u^* + 2 (L^T x - w)^T u^* \\
&= -(L^T x - w)^T (R_d + B_d^T P B_d)^{-1} (L^T x - w) \\
&= -x^T K_d^T (R_d + B_d^T P B_d) K_d x + 2 w^T K_d x + \text{const}
\end{aligned}
$$

i.e.
$$
J(x, u^*) = x^T \Big( Q_d + A_d^T P A_d - K_d^T (R_d + B_d^T P B_d) K_d \Big) x - 2 \Big( s_x + A_d^T v - K_d^T w \Big)^T x + \text{const}
$$

Matching the quadratic and linear terms, the new cost-to-go is $x^T P x - 2 v^T x$ with
$$
\begin{aligned}
P &\leftarrow Q_d + A_d^T P A_d - K_d^T (R_d + B_d^T P B_d) K_d \\
v &\leftarrow A_d^T v + s_x - K_d^T w \\
&= A_d^T v + s_x - K_d^T (B_d^T v + s_u) \\
&= (A_d - B_d K_d)^T v + (s_x - K_d^T s_u) \\
&= A_{cl,d}^T v + \sum_{i=0}^{d-1} (A_i - B_i K_d)^T q_{j+i}
\end{aligned}
$$

**Joseph-like form.** Expanding $A_{cl,d}^T P A_{cl,d}$ and using $K_d^T (R_d + B_d^T P B_d) = N_d + A_d^T P B_d$ recovers the Joseph-like form given previously (in the MPC $d$-step recursion section):
$$
\begin{aligned}
&Q_d + K_d^T R_d K_d - N_d K_d - K_d^T N_d^T + A_{cl,d}^T P A_{cl,d} \\
&= Q_d + K_d^T R_d K_d - N_d K_d - K_d^T N_d^T + (A_d - B_d K_d)^T P (A_d - B_d K_d) \\
&= Q_d + K_d^T R_d K_d - N_d K_d - K_d^T N_d^T + A_d^T P A_d - A_d^T P B_d K_d - K_d^T B_d^T P A_d + K_d^T B_d^T P B_d K_d \\
&= Q_d + A_d^T P A_d + K_d^T (R_d + B_d^T P B_d) K_d - (N_d + A_d^T P B_d) K_d - K_d^T (N_d + A_d^T P B_d)^T \\
&= Q_d + A_d^T P A_d + K_d^T (R_d + B_d^T P B_d) K_d - K_d^T (R_d + B_d^T P B_d) K_d - K_d^T (R_d + B_d^T P B_d) K_d \\
&= Q_d + A_d^T P A_d - K_d^T (R_d + B_d^T P B_d) K_d
\end{aligned}
$$

> Note: the identity $K_d^T (R_d + B_d^T P B_d) = N_d + A_d^T P B_d$ comes from the gain definition $K_d = (R_d + B_d^T P B_d)^{-1} (B_d^T P A_d + N_d^T)$ and the fact that $R_d + B_d^T P B_d$ is symmetric.

**Held $r$.** Replacing each $q_{j+i}$ with the block average $\bar{q}$ gives:
$$
s_x = \sum_{i=0}^{d-1} A_i^T \bar{q} = \Big( \sum_{i=0}^{d-1} A_i \Big)^T \bar{q} = G_d^T \bar{q}, \qquad s_u = \sum_{i=0}^{d-1} B_i^T \bar{q} = \Big( \sum_{i=0}^{d-1} B_i \Big)^T \bar{q} = H_d^T \bar{q}
$$
Substituting these into the $v$ update ($s_x - K_d^T s_u$ form) and the optimal input, we recover the previous held-$r$ expressions:
$$
\begin{aligned}
v &\leftarrow A_{cl,d}^T v + (G_d - H_d K_d)^T \bar{q} \\
u^* &= -K_d x + (R_d + B_d^T P B_d)^{-1} (B_d^T v + H_d^T \bar{q})
\end{aligned}
$$

**Held $x$.** Setting $A_i = I$, $B_i = 0$ inside the stage-cost sums (the transition $A_d$, $B_d$ is untouched) gives $Q_d = dQ$, $R_d = dR$, $N_d = 0$, $G_d = dI$, $H_d = 0$. This recovers the simplified recursion, and the feedforward reduces to $(dR + B_d^T P B_d)^{-1} B_d^T v$: the within-block term $H_d^T \bar{q}$ is exactly what the held-$x$ assumption neglects.

**Check: $d = 1$.** With $A_1 = A$, $B_1 = B$, $A_0 = I$, $B_0 = 0$: $Q_1 = Q$, $N_1 = 0$, $R_1 = R$, $s_x = q_j$, $s_u = 0$. The block recursion reduces exactly to the fine MPC recursion.

$\square$

### $P$ under decimation

**In the coarse recursion**, $P$ is still the fine-rate cost-to-go (the stage cost counted at every $T_s$ step), just under the constraint that $u$ is held across each block. This is why it carries no $d$ subscript: it is the same quantity, computed at block boundaries. A gain computed via $P$ is still a valid control law for the system at sampling rate `Ts`, just under held-u/held-x assumptions.

**In preview control**, decimation doesn't affect $P$ at all. The fine-rate feedback acts at every step, so $P$ stays at the DARE solution, and decimation only changes how $v$ is accumulated.

**In MPC control,** each coarse step computes $K_d$ from $P$ using the $d$-step matrices. As a control law, $K_d$ is only suited to an input held across the block: under held-$x$, a purely coarse recursion gives essentially `dlqr` at `d*Ts`, so for large $d$ it would be badly mismatched if applied at `T_s`. On the other hand, $P$ stays meaningful: it is the *fine-rate cost-to-go* (in pure-coarse steady state, $d$ times the `dlqr` solution at `d*Ts`), just inflated according to the held-u assumption. 

In the implementation, the coarse gains are never applied: they only propagate $P$ and $v$ backward. The applied gain comes from the fine law $K_k = (R + B^T P B)^{-1} B^T P A$ with the final $P$, computed via fine steps. The held-$u$ and held-$x$ approximations then enter only through $P$ as a conservative estimate of the far-horizon cost, keeping the input free to move in the nearest horizon (where it matters most).

### Held-x

**What the held-$x$ assumption does.** Take a block with $d = 4$. The transition across the block is exact in both cases; the difference is which states the stage cost is evaluated at:

    fine step:        j        j+1       j+2       j+3     |   j+4
                                                           |
    exact:            x  ───>  x_1  ───>  x_2  ───>  x_3  ──┼─>  x^+ = A_d x + B_d u
    stage cost at:    x        x_1        x_2        x_3    |
                                                           |
    held-x:           x  ───>  x_1  ───>  x_2  ───>  x_3  ──┼─>  x^+ = A_d x + B_d u
    stage cost at:    x        x          x          x      |

where $x_i = A_i x + B_i u$. Under held-$x$, the stage cost counts the block-start state $d$ times (hence $dQ$, $dR$, $d\bar{q}$) rather than following the state through the block, i.e. a left-rectangle approximation of the total block cost.

**Why this is reasonable.** The $d$-step dynamics are untouched, so only the cost within each block is approximated. This error scales with how far the state moves within a block – which is small for typical $d$, and located in the far-horizon, where the prediction is least accurate, and fidelity matters least. In testing, the exact version added complexity without improving performance.

**Parallel with a $d T_s$ discretisation.** Held-$x$ makes each block exactly *what a system sampled at* $d T_s$ *sees*: the ZOH transition $(A_d, B_d)$, with the state and input cost applied only at the block boundaries. Each coarse step is therefore exactly the Riccati step of the $d T_s$ system with scaled costs, $(A_d, B_d, dQ, dR)$; $K_d$ is that system's gain, and the scaling only changes the size of $P$, keeping it in the same units as the fine steps it's combined with. It replaces $d$ fine steps with one evaluation – but retaining $P$ in fine-rate units, so that coarse and fine steps can be mixed freely throughout the recursion.