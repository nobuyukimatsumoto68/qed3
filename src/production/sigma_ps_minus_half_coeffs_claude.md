# Four-point coefficients for $\sigma = \mathrm{PS} - \tfrac12$ (c-number contact subtraction)

Nf2 $g^2{=}0.5$ L1, $a_t{=}0.2$, massless, 400 cfg, distillation $N_v{=}24$ (exact).
Codes: `two_meson_cnumber_claude.py` (corrected), `diag_effmass_claude.py` (old, wrong).

## 1. The operator: $-\tfrac12$ is a c-number, not a kernel modification

$\sigma(x) = \bar\psi S\psi(x) - \tfrac12$. The $-\tfrac12$ is the **contact VEV** of the local density. From the
overlap propagator $D_\text{ov}^{-1} = \tfrac12\,\delta(x,y) + (\text{propagation})$ and the loop/anticommutation
sign,
$$
\langle\bar\psi_A\psi_B\rangle = -(D_\text{ov}^{-1})_{BA}\ \xrightarrow{\text{contact}}\ -\tfrac12\delta_{AB}.
$$
The diagram trace $D_S=+\mathrm{Tr}[\Phi\,\tau(s,s)] = -\langle\bar\psi S\psi\rangle$ carries the opposite (closed-loop)
sign, so subtracting $c_0 = \tfrac12\mathrm{Tr}\,\Phi$ from the *diagram* $D_S$ is what enforces $\langle\sigma\rangle=0$;
verified $\langle D_S - c_0\rangle = 0$ to $10^{-7}$. The $\tfrac12$ magnitude and the net normal-ordering are correct.

**Key point.** $-\tfrac12$ is attached to a **self-contraction** (a bilinear contracting with itself = the tadpole).
It is NOT a modification of the propagator that links two *different* fields. Therefore:
- it cancels the standalone **tadpole** $D_S$ and the disconnected vacuum;
- the equal-time $\tau(s,s)$ **inside** the connected loops $D'_S,V_S,S_S$ links two *different* $\sigma$'s and is
  a genuine part of the connected loop -- it stays **RAW** (contact included).

The old code subtracted $\tau(s,s)\to\tau(s,s)-\tfrac12 I$ in *every* leg (a kernel/propagator modification). That is
correct only for the standalone $D_S$; it wrongly modifies $D'_S,V_S,S_S$.

## 2. Prescription (tadpole-only)

In the 10 diagrams:
$$
D_S \to D_S - c_0\ \ (c_0=\tfrac12\mathrm{Tr}\,\Phi)\quad\text{on standalone tadpole factors only};\qquad
D'_S,\,V_S,\,S_S,\,C_S,\,T_S:\ \text{RAW } \tau .
$$
Equivalently, expand $\langle(\tilde P-c_0)^2(t)(\tilde P-c_0)^2(0)\rangle$ in raw correlators of $\tilde P=\bar\psi S\psi$;
the shifts cancel the tadpole VEVs. Checked on the composite one-point:
$\langle\sigma^2\rangle = \langle(D_S-c_0)^2 + D'_S^\text{raw}\rangle = \langle D_S^2\rangle - c_0^2 + \langle D'_S\rangle$,
matching both routes.

## 3. Corrected diagram set

Combinatorial weights $\{4,2,4,4,2,1,4,1,1,1\}$ are **unchanged** (they count contraction patterns). What changes is
the operator content: with $\delta D_S = D_S - c_0 \approx 0$ (std $10^{-7}$), the tadpole-carrying diagrams
C, D, G, H, I, J $\to 0$; the survivors are those with no standalone tadpole, evaluated with RAW $\tau$:
$$
\boxed{\ \langle O(t)O(0)\rangle \;\approx\; 2\big[\,4(-S_S) + 2(-T_S) + 2\,C_S^2 + D'_S(s)D'_S(t)\,\big]\ }
$$
surviving set $\{A,B,E,F\}$ with weights $\{4,2,2,1\}$, raw $\tau$.

## 4. What changes numerically (vs old $\tau-\tfrac12 I$)

| diagram | old | new (c-number) | note |
|---|---|---|---|
| A $(-S_S)$ | $m{=}0.53$, $|C(4)|{=}4{\times}10^{-7}$ | $m{=}0.36$, $|C(4)|{=}4.5{\times}10^{-5}$ | raw $\tau$: 100$\times$ larger, single-$\sigma$ |
| B $(-T_S)$ | $0.69$ | $0.69$ | unchanged (no equal-time legs) |
| E $(C_S^2)$ | $0.72$ | $0.72$ | unchanged (no equal-time legs) |
| F $(D'_S D'_S)$ | noise $10^{-8}$ | $|C(4)|{=}4.5{\times}10^{-6}$ | now real: $\langle D'^\text{raw}_S\rangle=0.514$ |
| C,D,G,H,I,J | $\sim10^{-9}$ | $\sim10^{-9}$ | tadpoles, vanish either way |

Composite one-point: physical $\langle\sigma^2\rangle = \langle D'^\text{raw}_S\rangle = \mathbf{0.514}$ (old code's
$-0.0197$ was the artifact). The old code was suppressing A and killing F -- the two connected diagrams with buried
equal-time legs.

## 5. Invariance of the connected cumulant

The old$\to$new difference is a set of $-\tfrac12 I$ insertions that each collapse an equal-time leg into a
**lower-point / disconnected** structure. The connected four-point cumulant
$\langle\sigma\sigma\sigma\sigma\rangle_c = \langle OO\rangle - \tfrac12\langle O\rangle^2 - 2\langle C_S\rangle^2$
subtracts exactly those. Hence the connected cumulant is **prescription-invariant**:
- the $a_t m_{\sigma\sigma}^\text{conn}\approx1.0$ benchmark (`sigma_connected_gevp_benchmark_claude.md`) **stands**;
- the F-connected 0++ test result is likewise robust (both prescriptions give noise).

The correction matters for the **full (non-connected) four-point**, the **per-diagram** decomposition, and the
**one-point** $\langle\sigma^2\rangle$ -- i.e. for building the correct operator basis for the GEVP, not for the
connected physics already extracted.

## 6. 0++ test (disconnected same-timeslice)

$\langle\delta D'_S(t)\,\delta D'_S(0)\rangle_c$, $\delta D'_S = D'_S - \langle D'_S\rangle$ (the local $\sigma^2$
density = fermionic 0++ candidate): comes out **negative, $|C/\text{err}|\sim1$--$2$**, no decaying plateau near
$a_t m_{F^2}=0.616$. So at this ensemble/stats the fermionic $\sigma^2$ density does **not** resolve a 0++ ->
favours **outcome 2** (0++ essentially gluonic, weak off-diagonal mixing), not outcome 1. Caveat: negative +
barely-resolved could be low stats or a reflection-positivity subtlety from the non-Hermitian $D_\text{ov}$; revisit
with more cfg. Plot: `figs/two_meson_cnumber_Fconn_Nf2_gsq0.500000at0.200000_claude.png`.

## 7. Status / next

- Corrected diagram code: `two_meson_cnumber_claude.py` (`diags_pair_cnumber`).
- Connected benchmark unaffected; no need to redo it.
- Open: the F-connected 0++ sign/stats; the FS-channel analogue (Step 3 of `sigma_sigma_next_steps_plan_claude.md`);
  then GEVP (Lanczos, gluonic mixing).
