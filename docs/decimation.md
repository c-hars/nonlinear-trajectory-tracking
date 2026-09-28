## Decimated recursion: glossary

Terminology, useful for understanding the functions in the repo and the math that follows.

- **Fine step:** a prediction step of duration `Ts`.
- **Coarse step:** a prediction step covering more time than a fine step (duration > `Ts`); lower model fidelity, faster compute. A coarse step has duration `d*Ts`, where `d` is the decimation factor; `d` consecutive fine steps are grouped and approximated as one coarse step.
- **Block:** each discrete-time step is over `d*Ts` seconds; 'block' is shorthand for this chunk of time.
- **Coarse matrices:** `(A_d, B_d, Qy_d, R_d)` – the system matrices `(A, B, Q, R)` after applying a decimation factor. They represent the same continuous-time dynamics and OCP formulation `(Ac, Bc, Qy_c, R_c)`, but discretised at a different (lower-resolution) sampling rate. Also referred to as `d`-step matrices.
- **Fine recursion:** a backward recursion processing one fine step at a time.
- **Coarse recursion:** a backward recursion processing one coarse step at a time.
- **Tier:** a contiguous sequence of steps sharing the same decimation factor.
  - **Far-horizon tier:** the final tier – farthest from the current time, and uses the largest decimation factor. Absorbs whatever horizon remains after the near-horizon tiers.
  - **Near-horizon tier:** any tier other than the far-horizon tier. Each holds `StepsPerTier` steps; each step spans a block of `d` fine steps.
- **Remainder fine steps:** fine steps that don't fit cleanly into the tier schedule; the prediction horizon of `M` steps might have a remainder. These are not discarded, instead they're placed nearest to the current time and processed by the fine recursion.
- **Preview grid:** the full arrangement across the preview horizon. Forward in time: remainder fine steps → near-horizon tiers → far-horizon tier → preview endpoint. The backward recursion processes this in reverse.

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
v \leftarrow F_d v + \Phi_d \bar{q}
$$

where
$F_d = F^d$ and $ \Phi_d = \sum_{i=0}^{d-1} F^i $.

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
for $i = 0, \dots, d-1$ and some $j \geq 0$. The exact $d$-step formula is then:

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

Then with $r$ held, the costate update becomes:
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
P &\leftarrow Q_d + K_d^T (R_d) K_d + A_{cl,d}^T P A_{cl,d} \\
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

---
