# Plan: Connected and Disconnected Current Correlator Measurements

## Physics / Goal Summary

This plan specifies new measurement programs that compute current-current two-point
functions using the exactly conserved overlap currents, defined on lattice links
$(n,n')$ at time $t$ (spatial links) and $(n,t)\to(n,t+1)$ (temporal links).

We compute **all** current types, not just the vector current:

- **Vector** $J_{V,\pm}$ (Eqs. 3.41--3.44): connected + disconnected.
- **Flavor** $J_{r,\pm}$ (Eqs. 3.49--3.52): connected only, equal to the vector connected
  piece.  Note the flavor currents **do not exist for $N_f = 2$** (the flavor group is
  trivial there); they are relevant only for $N_f = 4, 6$.
- **Axial** $J_{A,\pm}$ (Eqs. 3.45--3.48, both orderings): connected only.

We will be working on $N_f = 2, 4, 6$.  Since the connected piece is shared
(vector-connected $=$ flavor), the genuinely distinct computations are: the disconnected
trace (vector only), the connected piece (vector/flavor), and the axial pieces.

For the vector current $J_{V,\pm}$ the correlator has **both** a connected diagram
($C_{V,c}$) and a disconnected diagram ($C_{V,d}$); see Eqs. (3.41)--(3.44).  The two
diagrams are computed by **different programs with different workflows**:

- **Disconnected ($C_{V,d}$):** it factorizes,
  $C_{V,d}^{nn'}(0\to t) = T^{nn'}(0)\,T^{nn'}(t)$ with
  $T^{nn'}(t) \equiv \mathrm{tr}[D_\text{ov}^{-1}K^{nn'}(t)]$.  The measurement program
  (`jj_disc_claude.cu`) only **dumps the per-config, per-link / per-site traces**
  $T^{nn'}(t)$ and $T^{t,t+1}(n)$ -- the raw building blocks.  It does **not** form the
  correlator, the gauge average, the geometric weighting, or the spherical-harmonic
  projection.  Those are all done afterwards in a **separate post-processing script**.

- **Connected ($C_{V,c}$):** it does not factorize, so it is computed as a correlator
  in-program from a **single source** via a source-side vector $\phi=K^{(\text{src})}D_\text{ov}^{-1}\eta$
  and a sink-side vector $\psi=D_\text{ov}^{-\dagger}K^{(\text{snk})\dagger}\eta$, with
  $C_{V,c}=\langle\psi^\dagger\phi\rangle$ (there is no $\gamma_5$-hermiticity in this
  system, so the naive single-inversion one-end product is wrong -- see "Connected
  estimator" below).  Split into three programs by projection type
  (`_spproj` / `_tpproj` / `_ylmproj`); a `--current vector|axial` CLI flag picks which
  current each run computes (no inversion reuse across currents -- simplicity first).  The
  vector piece is identical to the flavor-current correlator (Eqs. (3.49)--(3.52)).

The final physics observables are the spatially projected correlator (Eq. (4.29)), the
temporally projected correlator (Eq. (4.32)), and the spherical-harmonic-projected
correlator (Eq. (4.36)).  The **geometric weight factors** ($A/\kappa^{(0)2}$,
$A\,Y_{\ell m}/\kappa^{(0)}$) are applied **in the compiled calculation**: the connected
programs include them directly and output the projected correlators; the disconnected dump
is raw traces, and the disc post-processing programs apply the weights.  Only the gauge
average $\langle\cdot\rangle_U$, vacuum subtraction (disc), the vector combination
$-C_{V,c}+C_{V,d}$, and fitting/plotting are left to the `.ipynb` notebooks.

The derivation and notation follow Sections 3.3 and 4 of
`/mnt/barracuda22/qed3/qed3int_v2.pdf` (Katz-Matsumoto), Eqs. (3.39)--(3.57) and
(4.27)--(4.36).

---

## Physics Background

### Notation (two-component, $N_f = 2$, $\xi$/$\eta$ basis)

The Wick contraction rules (Eqs. (3.39)--(3.40)), written with the fermionic average
$\langle\,\cdot\,\rangle_F$ as in the PDF, are

$$\langle \xi_{f,x_1}\,\eta_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-1})_{x_1 x_2}, \qquad
  \langle \eta_{f,x_1}\,\xi_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-\dagger})_{x_1 x_2}.$$

For a spatial link $(n,n')$ on timeslice $t$, the kernel is $K^{nn'} \equiv K^{wz}$
(two-component, Eq. (3.32) and its adjoint), implemented in
`includes/conserved_current_claude.h::ConservedCurrent::apply_k`.  The same `apply_k`
handles the temporal link $K^{t,t+1}(n)$ via the `std::pair<int,Idx>` overload.

The (+) branch of the vector current (two-component, Eq. (3.36)) is

$$J_{V,+}^{nn'}(t) = \eta_f^\dagger K^{nn'}(t) \xi_f.$$

### Vector correlator: connected + disconnected (Eqs. (3.41)--(3.44))

$$
\frac{2}{N_f}\langle J_{V,+}^{nn'}(0)\,J_{V,+}^{nn'}(t)\rangle_F
  = -\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]
  + \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\right]\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(t)\right]
$$

$$\equiv -C_{V,c}^{nn'}(0\to t;U) + C_{V,d}^{nn'}(0\to t;U). \tag{3.41, 3.42}$$

The connected piece (does **not** factorize):

$$C_{V,c}^{nn'}(0\to t;U) = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]. \tag{3.41}$$

The disconnected piece (**factorizes** into two single-time traces):

$$
C_{V,d}^{nn'}(0\to t;U)
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\right]
    \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(t)\right]
  = T^{nn'}(0)\,T^{nn'}(t),
$$

with the single-time trace

$$T^{nn'}(t) \equiv \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(t)\right]. \tag{3.42}$$

Note $T^{nn'}(t) = -J_{V,+}^{nn'}(t)/N_f$ from the formula
$J_{V,+} = -N_f\,\mathrm{tr}(K D_\text{ov}^{-1})$ already verified in
`check_conserved_current_claude.cu` Chunks 5b and 6.  The disconnected program dumps
exactly $T^{nn'}(t)$ (equivalently $L^{nn'}(t) = -N_f T^{nn'}(t)$); the choice of the
overall $-N_f$ factor is a post-processing convention.

### Symmetry under time reversal (Eqs. (3.43)--(3.44))

$$
\frac{2}{N_f}\langle J_{V,-}^{nn'}(0)\,J_{V,-}^{nn'}(t)\rangle_F
  = -\bigl(C_{V,c}^{nn'}(0\to t;U)\bigr)^* + \bigl(C_{V,d}^{nn'}(0\to t;U)\bigr)^*. \tag{3.43, 3.44}
$$

Consistency check: $C_{V,c}(t)^* = C_{V,c}(N_t - t)$ etc.

### Flavor current correlator (Eqs. (3.49)--(3.52))

For the flavor $SU(N_f/2)^2$ current (Eq. (3.38), (+) branch):

$$
\langle J_{r,+}^{nn'}(0)\,J_{r',+}^{nn'}(t)\rangle_F
  = -\delta_{rr'}\,\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]
  = -\delta_{rr'}\,C_{V,c}^{nn'}(0\to t;U). \tag{3.49, 3.50}
$$

Identical to the connected vector piece; `jj_conn_claude.cu` serves both.

### Axial correlator $C_A^{nn'}$ (Eqs. (3.45)--(3.48)) -- in scope

Both orderings are computed.  Forward:

$$
\frac{2}{N_f}\langle J_{A,+}^{nn'}(0)\,J_{A,-}^{nn'}(t)\rangle_F
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)(1-D_\text{ov}^\dagger)\,
      D_\text{ov}^{-\dagger}K^{nn'\dagger}(t)(1-D_\text{ov})\right]
  \equiv C_A^{nn'}(0\to t;U). \tag{3.45, 3.46}
$$

Backward (the $D\leftrightarrow D^\dagger$, $K\leftrightarrow K^\dagger$ swap):

$$
\frac{2}{N_f}\langle J_{A,-}^{nn'}(0)\,J_{A,+}^{nn'}(t)\rangle_F
  = \mathrm{tr}\!\left[D_\text{ov}^{-\dagger}K^{nn'\dagger}(0)(1-D_\text{ov})\,
      D_\text{ov}^{-1}K^{nn'}(t)(1-D_\text{ov}^\dagger)\right]
  = C_A^{nn'}(t\to 0;U). \tag{3.47, 3.48}
$$

Purely connected (no disconnected term); one leg $D_\text{ov}^{-1}$, the other
$D_\text{ov}^{-\dagger}$, with the Ginsparg-Wilson factors $(1-D_\text{ov})$,
$(1-D_\text{ov}^\dagger)$.  The Sec. 4 projections (4.29/4.32/4.36) apply to the axial
channel as well (the estimators are generic over current type).

### Stochastic estimators

$Z_2 \times Z_2$ noise throughout.

**Disconnected single-time trace** at link $(n,n')$, timeslice $t$ (with time-spin
dilution as in `disc_claude.cu`):

$$
T^{nn'}(t) = \mathrm{tr}(K^{nn'}(t)\,D_\text{ov}^{-1})
  \approx \frac{1}{N_\text{hits}} \sum_h \eta_h^\dagger K^{nn'}(t)\,\phi_h,
  \qquad D_\text{ov}\phi_h = \eta_h.
$$

This is the estimator already working in `check_conserved_current_claude.cu` Chunk 5b.
The program dumps $T^{nn'}(t)$ per config; the disconnected correlator
$\langle T(0)\,T(t)\rangle_U$ (gauge average, with any vacuum subtraction) is formed in
post-processing.

**Connected piece -- single source, source-side and sink-side vectors.**

There is **no $\gamma_5$-hermiticity** in this system: $D_\text{ov}^{-\dagger}\neq
\gamma_5 D_\text{ov}^{-1}\gamma_5$.  Consequently the naive single-inversion one-end product
fails.  With one solve $\phi=D_\text{ov}^{-1}\eta$ and $c_0(n)=\eta^\dagger K(n)\phi$:

$$
E\!\left[c_0(n_1)\,c_0(n_2)^*\right]
  = \mathrm{tr}[K(n_1)D_\text{ov}^{-1}]\,\bigl(\mathrm{tr}[K(n_2)D_\text{ov}^{-1}]\bigr)^*
  + \mathrm{tr}\!\left[K(n_1)D_\text{ov}^{-1}\,D_\text{ov}^{-\dagger}K(n_2)^\dagger\right],
$$

whose connected term carries $D_\text{ov}^{-1}D_\text{ov}^{-\dagger}$, not the
$D_\text{ov}^{-1}\!\cdots D_\text{ov}^{-1}$ of the target Eq. (3.41).

The correct estimator (single source $\eta$) builds a **source-side** vector with the
current applied *after* one inversion, and a **sink-side** vector with the adjoint current
applied *before* an adjoint inversion:

$$
\phi = K^{(\text{src})}\,D_\text{ov}^{-1}\eta,
\qquad
\psi = D_\text{ov}^{-\dagger}\,K^{(\text{snk})\dagger}\,\eta,
\qquad
\mathcal{E} = \psi^\dagger\phi .
$$

Since $\psi^\dagger = \eta^\dagger K^{(\text{snk})}D_\text{ov}^{-1}$,

$$
E[\mathcal{E}]
  = E\!\left[\eta^\dagger K^{(\text{snk})}D_\text{ov}^{-1}K^{(\text{src})}D_\text{ov}^{-1}\eta\right]
  = \mathrm{tr}\!\left[K^{(\text{snk})}D_\text{ov}^{-1}K^{(\text{src})}D_\text{ov}^{-1}\right]
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{(\text{snk})}D_\text{ov}^{-1}K^{(\text{src})}\right]
  = C_{V,c}.
$$

A single bilinear $\eta^\dagger M\eta$ with $M$ fixed -> $E[\mathcal{E}]$ is **exactly** the
trace; both legs are forward $D_\text{ov}^{-1}$, no contamination, no spurious $D^{-\dagger}$.
($K^\dagger$ is provided by `ConservedCurrent::apply_k_dag` (`conserved_current_claude.h:157`,
spatial wrapper `apply_k_dag_wz` `:230`) -- implemented, under verification.)

**Key efficiency -- source side free; one composite sink solve per sink time.**  Solve
$\phi_0 \equiv D_\text{ov}^{-1}\eta$ **once** per hit; the source current $K^{(\text{src})}$
is applied to $\phi_0$ *after* the inversion, so **every source insertion point and time is
free**.  The sink side inverts $D_\text{ov}^\dagger$ on a **composite** source (the
insertion-point sum folded in -- $Y_{\ell m}$ weights for ylm, $\zeta$-noise for the
diagonal; see "Projection dependence"), so it is **one $D^\dagger$ solve per sink time**, not
per insertion point.

**Insertion-point averaging with the $Y_{\ell m}$ weights (Eq. 4.36).**  Fold the
projection sums into the source and sink vectors:

$$
\Phi^{(\ell_1 m_1)}(s_0) = \sum_{n_1}\frac{A_{n_1}Y_{\ell_1 m_1}(\hat n_1)}{\kappa^{(0)}_{t,t+1}(n_1)}\,
  K^{t,t+1}(n_1,s_0)\,\phi_0,
$$

$$
\Psi^{(\ell_2 m_2)}(s) = D_\text{ov}^{-\dagger}\sum_{n_2}\frac{A_{n_2}Y_{\ell_2 m_2}(\hat n_2)}{\kappa^{(0)}_{t,t+1}(n_2)}\,
  K^{t,t+1\,\dagger}(n_2,s)\,\eta,
$$

$$
G^{t}_{(\ell_1 m_1),(\ell_2 m_2)}(t) = \Big\langle\,\frac{1}{N_t}\sum_{s_0}
  \Psi^{(\ell_2 m_2)}(s_0+t)^\dagger\,\Phi^{(\ell_1 m_1)}(s_0)\,\Big\rangle_U .
$$

$\Phi$ for all source times $s_0$ and all source channels is free from $\phi_0$; $\Psi$
needs one $D^\dagger$ solve per sink channel per sink time used.

**Projection dependence (cost; why three programs).**  In all cases the sink solve is on a
single **composite** source built by folding the projection into the sum over insertion
points -- one $D^\dagger$ solve per sink time, *not* per insertion point.

- **`_ylmproj` (Eq. 4.36):** the composite folds the deterministic $Y_{\ell m}$ weight --
  $n_\text{ch}$ $D^\dagger$ solves per sink time.
- **`_spproj` / `_tpproj` (Eqs. 4.29/4.32, link-/site-diagonal):** the diagonal sum is
  realized with a **stochastic insertion-point noise** $\zeta_a$ ($a$ = link / site,
  $E[\zeta_a\zeta_b^*]=\delta_{ab}$), so the composite

$$
\Phi^\zeta(s_0)=\sum_a \zeta_a\sqrt{w_a}\,K^a(s_0)\,\phi_0
\ (\text{free}),
\qquad
\Psi^\zeta(s)=D_\text{ov}^{-\dagger}\sum_a \zeta_a\sqrt{w_a}\,K^{a\dagger}(s)\,\eta
\ (\text{one solve}),
$$

  gives $E_{\zeta,\eta}[\Psi^\zeta(s)^\dagger\Phi^\zeta(s_0)] = \sum_a w_a\,
  \mathrm{tr}[K^a(s)D_\text{ov}^{-1}K^a(s_0)D_\text{ov}^{-1}]$, the diagonal projection.
  Cost: one $D^\dagger$ solve per sink time per $\zeta$-hit -- **no** $n_\text{links}$ /
  $n_\text{sites}$ inversion scaling.  Trade: extra $\zeta$ variance (off-diagonal terms
  average to zero), controlled by the number of $\zeta$-hits.

Because $\Phi$ gives *all* source times for free, the number of sink times at which $\Psi$
is solved sets both the cost and the source-time averaging: solving $\Psi$ at $n_s$ sink
times covers all $\Delta t$ with $n_s$ averaging samples (minimum 1).

All connected programs use all-to-all $Z_2\times Z_2$ noise (no time dilution) and report
the per-config projected correlator for offline gauge averaging.

**Axial estimator -- same source-side/sink-side split (Q1: confirm ordering).**

The axial correlator is purely connected, with Ginsparg-Wilson factors $(1-D_\text{ov})$,
$(1-D_\text{ov}^\dagger)$.  Target (forward, Eq. 3.45):

$$
C_A(0\to t) = \mathrm{tr}\!\left[D_\text{ov}^{-1}K(0)(1-D_\text{ov}^\dagger)\,
  D_\text{ov}^{-\dagger}K^\dagger(t)(1-D_\text{ov})\right].
$$

Use the same split $\mathcal{E}_A = \psi_A^\dagger\,\phi_A$, choosing the cut between the two
inverse propagators so that the **source side reuses a base inversion** (kernel applied
after) and the **sink side** carries the inversion (kernel applied before):

$$
\phi_A = D_\text{ov}^{-\dagger}\,K^\dagger(t)\,(1-D_\text{ov})\,\eta
\quad(\text{sink side; one } D^\dagger \text{ solve}),
\qquad
\psi_A = (1-D_\text{ov})\,K^\dagger(0)\,\rho_0
\quad(\text{source side; free}),
$$

with the base inversion $\rho_0 \equiv D_\text{ov}^{-\dagger}\eta$.  Then $\psi_A^\dagger =
\eta^\dagger D_\text{ov}^{-1}K(0)(1-D_\text{ov}^\dagger)$, so

$$
E[\psi_A^\dagger\phi_A]
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K(0)(1-D_\text{ov}^\dagger)\,
      D_\text{ov}^{-\dagger}K^\dagger(t)(1-D_\text{ov})\right]
  = C_A(0\to t).
$$

Note both inversions appear as $D_\text{ov}^{-\dagger}$ in this organization (the left
$D_\text{ov}^{-1}$ of the trace becomes $D_\text{ov}^{-\dagger}$ inside $\psi_A^\dagger$).
So the axial **source side** uses the base solve $\rho_0=D_\text{ov}^{-\dagger}\eta$
(kernel + GW factor applied after -> free over insertion points), and the axial **sink
side** $\phi_A$ is a composite $D^\dagger$ solve per sink time (same $Y_{\ell m}$ / $\zeta$
folding as the vector) -- with $K^\dagger$ and a $(1-D_\text{ov})$ on the source $\eta$.
The backward $C_A(t\to 0)$ (Eq. 3.47) is the $K\leftrightarrow K^\dagger$, $D\leftrightarrow
D^\dagger$ mirror, with source side using $\phi_0=D_\text{ov}^{-1}\eta$.

**Q1 (resolved):** the axial formula above (Eqs. 3.45/3.47) is correct and is to be
implemented **as written** -- forward $C_A(0\to t)$ has $K^\dagger(t)$ on the sink side
(behind the inversion) and $K(0)$ on the source side, with $(1-D_\text{ov}^\dagger)$
adjacent to $K(0)$ and $(1-D_\text{ov})$ adjacent to $K^\dagger(t)$ (mirror for backward).

**Current selection by CLI (Q2: decided -- simplicity over reuse).**

We do **not** try to reuse inversions across current types.  Each connected program run
computes **one** current type, selected by a CLI flag `--current vector|axial`:

- `--current vector`: the vector-connected piece $C_{V,c}$ (also the flavor-connected piece).
- `--current axial`: the two axial correlators $C_A(0\to t)$ and $C_A(t\to 0)$.

The vector and axial estimators each follow the source-side/sink-side construction above
(vector: $\phi=K^{(\text{src})}D_\text{ov}^{-1}\eta$, $\psi=D_\text{ov}^{-\dagger}K^{(\text{snk})\dagger}\eta$;
axial: the $\phi_A,\psi_A$ split with the GW factors), computed straightforwardly within
the selected branch.  No shared base solves between currents -- prioritize simple,
readable code.

### Section 4 estimators (geometric weights applied in the compiled calculation)

These are the physics observables.  The geometric weights below are applied **in the
compiled programs**: the connected programs fold them directly into the projected
correlator; the disc post-processing programs apply them to the raw trace dump.  The
notebooks then only gauge-average and fit.

Spatially projected, equal-link spatial average (Eq. (4.29)):

$$
G_s(t) \approx \frac{1}{4\pi}\sum_{\langle n,n'\rangle}
  \frac{A_{nn'}}{\kappa^{(0)\,2}_{nn'}}\,
  \langle J^{nn'}(0)\,J^{nn'}(t)\rangle .
$$

Temporally projected, site-diagonal (Eq. (4.32)):

$$
G_t(t) \approx \frac{1}{4\pi}\sum_{n}
  \frac{A_n}{\kappa^{(0)\,2}_{t,t+1}(n)}\,
  \langle J^{0,1}(n)\,J^{t,t+1}(n)\rangle .
$$

Spherical-harmonic-projected temporal correlator (Eq. (4.36)) -- needs the full
site$\times$site temporal correlator:

$$
G^{t}_{\ell_1 m_1;\ell_2 m_2}(t) \approx \frac{1}{4\pi}\sum_{n_1,n_2}
  \frac{A_{n_1}\,Y_{\ell_1 m_1}(\hat n_1)}{\kappa^{(0)}_{t,t+1}(n_1)}\,
  \frac{A_{n_2}\,Y_{\ell_2 m_2}(\hat n_2)}{\kappa^{(0)}_{t,t+1}(n_2)}\,
  \langle J^{0,1}(n_1)\,J^{t,t+1}(n_2)\rangle .
$$

In all three, $\langle J(0)J(t)\rangle = -C_{V,c} + C_{V,d}$ for the vector current
(Eqs. (3.41)--(3.42)).  The conformal-tower prediction (Eqs. (4.34)--(4.35)) lives in
the $Y_{\ell m}$ channels of Eq. (4.36): the $\ell=0$ channel vanishes (charge
conservation), the $\ell=1$ channel decays as $e^{-2t}$, the $\ell=2$ channel as
$e^{-3t}$.

### Geometric quantities (already in the code)

All weights needed by the Section 4 estimators are precomputed by the lattice classes:

| Symbol | Meaning | Code accessor |
|---|---|---|
| $\kappa^{(0)}_{nn'}$ | free-field link coupling | `DiracOp::kappa[il]` (`dirac_simp.h:355-360`, $= 2\,\texttt{link\_volume}/\texttt{ell}/\texttt{mean\_ell}$) |
| $A_{nn'}$ | dual area of edge $(n,n')$ (both triangles) | `lattice.link_volume[il]` (`s2n_simp.h:29`) |
| $\ell_{nn'}$ | primal link length | `lattice.ell[il]` (`s2n_simp.h:17`) |
| $A_n$ | dual (Voronoi) area of site $n$ | `lattice.dual_areas[s]` (`s2n_simp.h:38`) |
| $Y_{\ell m}(\hat n)$ | (real) spherical harmonic at site $n$ | `lattice.GetYlm(s,l,m)` after `lattice.UpdateYlm(l_max)` (`s2.h:55-56`) |

Eq. (4.36) uses **real** spherical harmonics $Y_{\ell,m}$ without the phase factor (see
the PDF Appendix C); `GetYlm` returns the complex $Y_\ell^m$, so the real combination is
formed where the weights are built -- i.e. in the compiled programs (connected and disc
post-proc), not in the notebooks.

### Lattice-to-continuum normalization (Eqs. (3.54)--(3.57))

The conserved currents do not renormalize, so the free-field couplings $\kappa^{(0)}$
relate the lattice operators to the continuum (Eq. (3.54)).  The triangle identity
(Eq. (3.56)) gives the spatial face projection (Eq. (3.57)) that underlies Eq. (4.29).
This geometric prefactoring is applied in the compiled programs (connected directly, disc
via the post-proc programs); the notebooks receive already-projected correlators.

---

## References

- A. Katz, N. Matsumoto, "Notes on QED3 (interacting theory)",
  `/mnt/barracuda22/qed3/qed3int_v2.pdf` --- Sec. 3.3 (lattice current correlators,
  Eqs. 3.39--3.57), Sec. 3.2.2 (two-component $K^{wz}$, Eq. 3.32), Sec. 4.2--4.3 (CFT
  estimators, Eqs. 4.27--4.36), Appendix C (real spherical harmonics).
  Primary physics reference.

---

## Scope and Constraints

**New files to create:**

| File | Role |
|---|---|
| `src/both_3d/jj_disc_claude.cu` | **Single** dump program: per-config kernel traces for *all* kernels -- $T^{nn'}(t)$ (spatial links) and $T^{t,t+1}(n)$ (temporal sites).  No projection; one file. |
| `src/both_3d/jj_conn_spproj_claude.cu` | Connected, spatial projection (Eq. 4.29): spatial-link kernels |
| `src/both_3d/jj_conn_tpproj_claude.cu` | Connected, temporal site-diagonal projection (Eq. 4.32): temporal-link kernels |
| `src/both_3d/jj_conn_ylmproj_claude.cu` | Connected, temporal $Y_{\ell m}$ projection (Eq. 4.36): temporal-link kernels, factorized weighted sum |
| `src/both_3d/jj_disc_spproj_claude.cu` | Disc post-processing, spatial projection (Eq. 4.29) |
| `src/both_3d/jj_disc_tpproj_claude.cu` | Disc post-processing, temporal site-diagonal projection (Eq. 4.32) |
| `src/both_3d/jj_disc_ylmproj_claude.cu` | Disc post-processing, temporal $Y_{\ell m}$ projection (Eq. 4.36) |

**Program-structure (decided):** the disconnected side has **one** dump program
`jj_disc_claude.cu` (the trace $T$ is the shared building block) plus **three**
post-processing programs `jj_disc_*proj_claude.cu` that read the dump and produce the three
projected correlator types in the same per-config layout as the connected programs; the
connected side is **three** measurement programs `jj_conn_*proj_claude.cu`.  All
projected outputs are then gauge-averaged / fit / plotted in the `.ipynb` notebooks.

**Existing files used (not modified):**

| File | Role |
|---|---|
| `includes/conserved_current_claude.h` | `ConservedCurrent::apply_k` (kernel $K^{wz}$, verified); `apply_k_dag` for $K^\dagger$ (`:157`, spatial wrapper `apply_k_dag_wz` `:230`) -- implemented, under verification |
| `includes/valence_claude.h` | `FermionVector`, `fill_z2_source`, `dag`, `time_spin_dilution` |
| `includes/overlap.h` | `Overlap<WilsonDirac>`, `adj_deviceAsyncLaunch`, `DHD_deviceAsyncLaunch` |
| `includes/s2n_simp.h` | `link_volume`, `ell`, `dual_areas`, `mean_ell` (geometric weights) |
| `qfe_mod/include/s2.h` | `UpdateYlm`, `GetYlm` (spherical harmonics) |
| `disc_claude.cu` | Template: dilution, CG solve chain, HDF5 output, resume logic |
| `check_conserved_current_claude.cu` | Template: `apply_k` + `dag` inner-product pattern |

**Compile-time parameters (inherit from `disc_claude.cu`):**

```cpp
namespace Comp {
  constexpr int N_REFINE = 1;
  constexpr int Nt       = 128;
  constexpr int NPARALLEL_DUPDATE = 1;
  const double TOL_INNER = 1.0e-15;
  const double TOL_OUTER = 1.0e-14;
}
```

**Runtime parameters (getopt_long, mirroring `disc_claude.cu`):**

| Flag | Default | Meaning |
|---|---|---|
| `--gsq` | 8.0 | Wilson coupling $g^2$ |
| `--Nf` | 2 | Number of sea quark flavors |
| `--nu0` | 1.0 | Sea quark asymmetry $\nu_0$ |
| `--nu1` | 1.0 | Valence quark asymmetry $\nu_1$ |
| `--nhits` | 1 | Number of $Z_2\times Z_2$ noise hits per config (production: 1, same as `disc_claude.cu`) |
| `--t_block` | 8 | Timeslices per dilution block (production: 8, same as `disc_claude.cu`; disc dump only) |
| `--current` | `vector` | Connected programs only: `vector` or `axial` -- which current to compute |

**Production settings:** `nhits=1`, `t_block=8` -- identical to the existing
`disc_claude.cu` (`disc_claude.cu:164-165`), with time-spin dilution
(`interval = Nt/t_block` blocks $\times$ `NS` spins per hit).

**Out of scope:**

- The post-processing script body (Eqs. 4.29/4.32/4.36 assembly) is listed as a separate
  artifact; its detailed plan is deferred until the dump format is fixed and validated.
- $L=2$ (`N_REFINE=2`) variant: same code, different `N_REFINE`.

---

## Data Structures

### Reuse: `FermionVector` (`includes/valence_claude.h`)

Host-pinned `field[N]` array of `Complex`.  Methods used:

- `fill_z2_source(rng)` --- all-to-all $Z_2\times Z_2$ noise (connected program)
- `time_spin_dilution(rng, t_s, t_block, spin)` --- diluted source (disconnected program,
  mirroring `disc_claude.cu:293`)
- `dag(b)` --- host inner product $\sum_i a_i^* b_i$

### Reuse: `ConservedCurrent<Fermion>` (`includes/conserved_current_claude.h`)

```cpp
// TEMPLATE: check_conserved_current_claude.cu:879
kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s, lk});       // spatial link  K
kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,Idx>{n, n_tnext});       // temporal link K
// adjoint kernel K^dag (same signature; conserved_current_claude.h:157 / :230)
kop.apply_k_dag(d_kdag_eta, d_eta_dev, U, std::pair<int,Idx>{n, n_tnext});
CUDA_CHECK(cudaDeviceSynchronize());
CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
const Complex c = eta.dag(fv_Kphi);
```

### New: per-config trace / correlator buffers

`jj_disc_claude.cu` (raw single-time traces, **no** correlator formed):

```cpp
// spatial: one complex per (spatial link, timeslice)
std::vector<std::vector<Complex>> T_spatial(base.n_links, std::vector<Complex>(Nt, 0.0));
// temporal: one complex per (site, timeslice)
std::vector<std::vector<Complex>> T_temporal(base.n_sites, std::vector<Complex>(Nt, 0.0));
```

Connected programs (source-side / sink-side vectors; correlator formed in-program,
per-config result dumped for offline gauge averaging):

```cpp
// base inversions (reused per hit; source sides apply kernels after these)
FermionVector phi0;   // = Dov^{-1} eta      (vector source side; axial backward source side)
FermionVector rho0;   // = Dov^{-dag} eta    (axial forward source side)
// source-side vectors (free over insertion points): Phi(s0) = sum_n w_n K(n,s0) phi0
// sink-side vector (one Dov^{-dag} solve per sink time, composite): Psi(s) = Dov^{-dag} sum_n w_n Kdag(n,s) eta
FermionVector Phi, Psi, Ksrc_phi, Kdag_eta;
// per-config projected correlators, per Delta-t.  sp/tp are already summed over
// links/sites (composite + zeta noise), so they are single length-Nt arrays; ylm keeps
// the (l,m) channel-pair matrix.
std::vector<Complex>              Cc_sp (Nt, 0.0);   // G_sp(Delta t)
std::vector<Complex>              Cc_tp (Nt, 0.0);   // G_tp(Delta t)
std::vector<std::vector<Complex>> Cc_ylm(n_ch*n_ch, std::vector<Complex>(Nt, 0.0));
FermionVector zeta_buf;  // insertion-point noise scratch (sp/tp)
```

---

## Ordered Implementation Chunks

### Chunk 1 --- `jj_disc_claude.cu`: boilerplate and construction

**What:** Copy boilerplate of `disc_claude.cu` (includes, `Comp` namespace, type aliases,
`ParseArgs`, GPU setup, lattice/DW/Overlap/Gauge construction).  Add
`#include "includes/conserved_current_claude.h"` and construct
`ConservedCurrent<Fermion> kop(D);`.

**Files:** `jj_disc_claude.cu`

```cpp
// TEMPLATE: disc_claude.cu:1-88   (includes, Comp namespace, type aliases)
// TEMPLATE: disc_claude.cu:156-229 (main: GPU setup, Base, DW, Fermion, Gauge)
// TEMPLATE: check_conserved_current_claude.cu:226-228 (ConservedCurrent construction)
```

---

### Chunk 2 --- `jj_disc_claude.cu`: enumerate links/sites and define dump layout

**What:** Build two enumerations:

```cpp
// spatial links: (timeslice s, base link lk)
for(int s = 0; s < Nt; s++)
  for(const auto& lk : base.links) spatial_list.push_back({s, lk});
// temporal links: (site n, timeslice s) -> link (n,s)->(n,s+1)
for(int s = 0; s < Nt; s++)
  for(int n = 0; n < base.n_sites; n++) temporal_list.push_back({n, s});
```

Output directory and HDF5 path (mirror `disc_claude.cu:206-214`):

```
jj_disc_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/jj_disc_trace.{k}.h5
```

**Dump layout (raw traces, no correlator):**

- `"spatial/{ix0}/{ix1}/real"`, `".../imag"` --- length `Nt` per spatial base link.
- `"temporal/{n}/real"`, `".../imag"` --- length `Nt` per site.

Resume sentinel: check existence of the last expected dataset before processing each $k$
(mirror `disc_claude.cu:272-278`).

**Files:** `jj_disc_claude.cu`

---

### Chunk 3 --- `jj_disc_claude.cu`: per-config stochastic loop (dump traces)

**What:** For each gauge config $k$, for each hit $h$, for each dilution block
$(t_s, \text{spin})$:

1. `eta.time_spin_dilution(rng, t_s, t_block, spin)`.
2. Solve $D\phi = \eta$ via normal equations:

```cpp
// TEMPLATE: disc_claude.cu:294-295
op_DH.from_cpu<N>(DH_eta.field, eta.field);
op_DHD.solve<N>(phi.field, DH_eta.field);
CUDA_CHECK(cudaMemcpy(d_phi_dev, reinterpret_cast<const CuC*>(phi.field), N*CD, H2D));
```

3. For each spatial link $(s, lk)$: apply $K^{nn'}$, accumulate
   $T_\text{spatial}[lk][s] \mathrel{+}= \eta^\dagger K^{nn'}(s)\phi$.
4. For each temporal link $(n, s)$: apply $K^{t,t+1}$, accumulate
   $T_\text{temporal}[n][s] \mathrel{+}= \eta^\dagger K^{t,t+1}(n)\phi$.

The dilution mask zeroes the contribution outside the active timeslices/spin, so summing
over all blocks reconstructs the full trace at each $(link, s)$.

5. After all blocks/hits, divide by `nhits`.  **Do not** form the correlator
   $\langle T(0)T(t)\rangle$ here -- dump $T_\text{spatial}$ and $T_\text{temporal}$ as-is.

**Files:** `jj_disc_claude.cu`

```cpp
// TEMPLATE: disc_claude.cu:288-299 (hit loop, dilution, solve, accumulate)
// TEMPLATE: check_conserved_current_claude.cu:870-885 (apply_k + dag inner product)
```

**Efficiency note:** each `apply_k` launches `D.size-1` CG solves (one per Zolotarev
pole).  Per solve there are `Nt*base.n_links` spatial + `Nt*base.n_sites` temporal
`apply_k` calls.  Sequential first; profile before optimizing (see Efficiency Decisions).

---

### Chunk 4 --- `jj_disc_claude.cu`: HDF5 dump and cleanup

**What:** Write the raw trace arrays to HDF5 with the keys of Chunk 2
(`HighFive::File::Truncate`, one file per $k$).  Cleanup `cudaFree` / `deallocate`.

```cpp
// TEMPLATE: disc_claude.cu:303-309 (createDataSet real/imag), disc_claude.cu:285-286 (Truncate)
```

**Files:** `jj_disc_claude.cu`

---

### Chunk 5 --- connected programs: boilerplate and construction

**What:** Same boilerplate as Chunk 1, for the three connected programs
`jj_conn_spproj_claude.cu`, `jj_conn_tpproj_claude.cu`, `jj_conn_ylmproj_claude.cu`.
All-to-all noise (no `--t_block`).  Add the `--current vector|axial` flag (selects which
current the run computes; a simple `if`/branch in the loop, no inversion reuse across
currents).  The three programs share the construction and the connected estimator (single
source; source-side vector free from $\phi_0=D_\text{ov}^{-1}\eta$, sink-side vector
$\psi=D_\text{ov}^{-\dagger}K^{(\text{snk})\dagger}\eta$ costs the solves; see "Connected
piece" in Physics Background).  They differ only in the projection and the output shape:

| Program | Links | Eq. | Sink-side composite solve (one per sink time) | Output |
|---|---|---|---|---|
| `_spproj` | spatial $(n,n')$ | 4.29 | $\zeta$-noise: $D_\text{ov}^{-\dagger}\sum_l \zeta_l\sqrt{w_l}\,K^{l\dagger}(s)\eta$ | $G_{sp}(t)$, length $N_t$ |
| `_tpproj` | temporal, site-diagonal | 4.32 | $\zeta$-noise: $D_\text{ov}^{-\dagger}\sum_n \zeta_n\sqrt{w_n}\,K^{t,t+1\dagger}(n,s)\eta$ | $G_{tp}(t)$, length $N_t$ |
| `_ylmproj` | temporal, off-diagonal | 4.36 | per channel: $D_\text{ov}^{-\dagger}\sum_n w\,K^{t,t+1\dagger}(n,s)\eta$ | $G^t_{\ell_1 m_1;\ell_2 m_2}(t)$, matrix |

**Files:** the three `jj_conn_*proj_claude.cu`

```cpp
// TEMPLATE: disc_claude.cu:1-229, check_conserved_current_claude.cu:226-228
```

---

### Chunk 6 --- connected programs: per-config stochastic loop

**What:** For each gauge config $k$, for each hit $h$ (single source $\eta$; see "Connected
piece" for $E[\psi^\dagger\phi]=C_{V,c}$).  Solve the **source-side base inversion once**,
apply source kernels after it (free), and pay $D^\dagger$ solves only on the sink side.

1. `eta.fill_z2_source(rng)` (all-to-all).
2. Base inversion: $\phi_0 = D_\text{ov}^{-1}\eta$ (normal equations, as Chunk 3).

`_tpproj` (Eq. 4.32, site-diagonal) -- composite over sites with insertion-point noise
$\zeta$, so the diagonal sum needs only one sink solve per sink time (no per-site solves):

```cpp
zeta = draw_z2(base.n_sites);                  // insertion-point noise, E[zeta_a zeta_b*]=delta
// source side (free): composite Phi(s0) = sum_n zeta_n sqrt(w_n) K(n,s0) phi0
for(s0){ Phi[s0]=0; for(n) axpy(Phi[s0], zeta[n]*sqrt(w_tp[n]), kop.apply_k(phi0,U,{n,n_tnext(n,s0)})); }
// sink side (one Ddag^{-1} solve per sink time s)
for(int s : sink_times){
  b = 0; for(n) axpy(b, zeta[n]*sqrt(w_tp[n]), kop.apply_k_dag(eta,U,{n,n_tnext(n,s)}));  // K^dag implemented
  op_DdagInv_solve(Psi, b);                     // Psi = Dov^{-dag} sum_n zeta_n sqrt(w_n) Kdag(n,s) eta
  for(s0) Cc_tp[ (s-s0+Nt)%Nt ] += Psi.dag(Phi[s0]);   // E_{zeta,eta}[.] = G_tp(Delta t)
}
```

`_spproj` (Eq. 4.29): identical with spatial-link kernels `std::pair<int,BaseLink>{s, lk}`
over `base.links` and $\zeta$ over links.  (Average over both $\eta$ and $\zeta$ hits.)

`_ylmproj` (Eq. 4.36, off-diagonal): fold the $Y_{\ell m}$ weights into the source-side
$\Phi$ (free) and the sink-side composite source (one solve per sink channel/time):

```cpp
// source side (free): Phi[ch1][s0] = sum_n w(ch1,n) K(n,s0) phi0
for(ch1,s0){ Phi[ch1][s0]=0; for(n) axpy(Phi[ch1][s0], w_ylm(ch1,n), kop.apply_k(phi0,U,{n,n_tnext(n,s0)})); }
// sink side: Psi[ch2][s] = Dov^{-dag} ( sum_n w(ch2,n) Kdag(n,s) eta )
for(int s : sink_times)
  for(int ch2=0; ch2<n_ch; ch2++){
    b = 0; for(n) axpy(b, w_ylm(ch2,n), kop.apply_k_dag(eta,U,{n,n_tnext(n,s)}));
    op_DdagInv_solve(Psi, b);                                // 1 solve per (ch2,s)
    for(ch1,s0) Cc_ylm[ch1*n_ch+ch2][ (s-s0+Nt)%Nt ] += Psi.dag(Phi[ch1][s0]);
  }
```

3. Average over $\eta$ hits (and, for `_spproj`/`_tpproj`, $\zeta$ hits); divide by the hit
   count.  The source-time average over $s_0$ is built in; each sink time $s$ in
   `sink_times` contributes all $\Delta t$, so $n_s$ sink times give $n_s$ source-time
   samples.
4. Write the per-config projected correlators (`Cc_sp` / `Cc_tp` / `Cc_ylm`) to HDF5; gauge
   averaging is offline.

`op_DdagInv_solve` is the $D_\text{ov}^{-\dagger}$ solve (analogous normal-equation chain to
Chunk 3, with $D_\text{ov}D_\text{ov}^\dagger$).

**Cost:** source side is free (one base solve + kernel applications).  Sink side is **one
$D^\dagger$ solve per sink time** for all three (a composite source): `_ylmproj` per channel,
`_spproj`/`_tpproj` per $\zeta$-hit.  No $n_\text{links}$ / $n_\text{sites}$ inversion
scaling -- the diagonal sum comes from the $\zeta$-noise average (Q2 resolution).

**Files:** the three `jj_conn_*proj_claude.cu`

```cpp
// TEMPLATE: disc_claude.cu:231-237,294-295 (CG solve chain for both inversions)
// TEMPLATE: check_conserved_current_claude.cu:870-885 (apply_k, dag)
```

---

### Chunk 7 --- connected programs: HDF5 output

**What:** Each program writes its per-config projected correlator to HDF5.  Output
directories (one per program):

```
jj_conn_spproj_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}/jj_conn_sp.{k}.h5
jj_conn_tpproj_..._/jj_conn_tp.{k}.h5
jj_conn_ylmproj_..._/jj_conn_ylm.{k}.h5
```

The filename also encodes `--current` (`vector` / `axial`).  Datasets: `_spproj` ->
$G_{sp}(t)$, length $N_t$ (already summed over links via the $\zeta$ composite); `_tpproj`
-> $G_{tp}(t)$, length $N_t$; `_ylmproj` -> per `(ch1,ch2)` channel pair, length $N_t$.
For `--current axial`, two such datasets ($0\to t$ and $t\to 0$).  Resume logic as in
Chunks 2/4.

**Files:** the three `jj_conn_*proj_claude.cu`

---

### Chunk 8 --- Validation

**What:** Compile-time-guarded self-checks (`if constexpr (N_REFINE==1 && Nt<=4)`, as in
`check_conserved_current_claude.cu:643`):

- `jj_disc`: compare stochastic $T^{nn'}(t)$ against the exact basis-vector trace.
- connected: verify $C_{V,c}(\Delta t=0) \geq 0$ (positivity at zero separation); optionally
  compare the $\langle\psi^\dagger\phi\rangle$ estimate against an exact
  $\mathrm{tr}[D^{-1}KD^{-1}K]$ at `Nt<=4`.

**Files:** guarded blocks in `jj_disc_claude.cu` / the `jj_conn_*proj_claude.cu` (no new file).

---

### Chunk 9 --- disc post-processing: setup, geometry, weights

**What:** Three `.cu`/`.cc` programs (one per projection type, mirroring the connected
split) that read the `jj_disc` trace dump and produce **per-config projected disconnected
correlators in the same layout as the connected programs** (Chunk 7).  Gauge averaging,
vacuum subtraction, fitting, and plotting are done downstream in the `.ipynb` notebooks.

| File | Eq. | Reads | Per-config output |
|---|---|---|---|
| `jj_disc_spproj_claude.cu` | 4.29 | spatial traces $T^{lk}(s)$ | per spatial link, length $N_t$ |
| `jj_disc_tpproj_claude.cu` | 4.32 | temporal traces $T^{t,t+1}(n,s)$ | per site, length $N_t$ |
| `jj_disc_ylmproj_claude.cu` | 4.36 | temporal traces | per channel-pair $(ch_1,ch_2)$, length $N_t$ |

(Compiled `.cu` if the lattice/geometry headers pull in CUDA; otherwise `.cc`.)

This chunk: reconstruct the lattice geometry (same construction as the measurement
programs -- read omega/alpha `.dat`, build `QfeLatticeS2`, `UpdateYlm(l_max)`), and build
the weight arrays:

$$
w^{sp}_{lk} = \frac{A_{lk}}{\kappa^{(0)\,2}_{lk}},\qquad
w^{tp}_{n} = \frac{A_{n}}{\kappa^{(0)\,2}_{t,t+1}(n)},\qquad
w^{ylm}_{n,(\ell m)} = \frac{A_{n}\,Y_{\ell m}(\hat n)}{\kappa^{(0)}_{t,t+1}(n)} .
$$

Enumerate the trace-dump files `jj_disc_trace.{k}.h5` over the config range; define the
output dir / HDF5 path mirroring Chunk 7.

**Files:** the three `jj_disc_*proj_claude.cu`

```cpp
// TEMPLATE: disc_claude.cu:156-229 (lattice/geometry construction), s2.h:758 (UpdateYlm)
```

---

### Chunk 10 --- disc post-processing: per-config projection

**What:** For each config $k$, read the dumped traces and form the per-config building
blocks (source-time average over $s_0$, wrapped mod $N_t$).  The disconnected correlator
factorizes ($C_{V,d}=T(0)\,T(t)$), so these are products of the same config's traces:

The **geometric weights are applied here, in-program** (built in Chunk 9); only the gauge
average / vacuum subtraction / fit are left to the notebooks.

`_spproj` / `_tpproj` (diagonal in link/site): apply $w^{sp}_{lk}$ (resp. $w^{tp}_n$) to the
per-config product, summed over links/sites:

$$
P^{sp}(t) = \sum_{lk} w^{sp}_{lk}\,\frac{1}{N_t}\sum_{s_0} T^{lk}(s_0)\,T^{lk}(s_0+t) .
$$

For correct vacuum subtraction the diagonal estimators must subtract $\langle T\rangle_U$
**per link/site before** the weighted sum, so the per-link/site single traces $T^{lk}(s)$
(weighted with $w^{sp}_{lk}$) are also output per config.

`_ylmproj` (factorized -- fold the $Y_{\ell m}$ weight first, then no site$\times$site
matrix):

$$
\tilde T^{(\ell m)}(s) = \sum_n w^{ylm}_{n,(\ell m)}\,T^{t,t+1}(n,s),
\qquad
D^{ylm}_{(\ell_1 m_1),(\ell_2 m_2)}(t) = \frac{1}{N_t}\sum_{s_0}
  \tilde T^{(\ell_1 m_1)}(s_0)\,\bigl(\tilde T^{(\ell_2 m_2)}(s_0+t)\bigr)^* .
$$

The (already-weighted) projected singles $\tilde T^{(\ell m)}(s)$ are also output per config
for the (factorized, per-channel) vacuum term.

**Files:** the three `jj_disc_*proj_claude.cu`

---

### Chunk 11 --- disc post-processing: output and notebook handoff

**What:** Write the per-config (already weighted) products and single-traces to HDF5 in the
connected-like layout (Chunk 7), one file per config $k$:

```
jj_disc_spproj_Nf{Nf}_gsq{gsq}.../jj_disc_sp.{k}.h5   # P^sp(t); per-link weighted single w*T^lk(s)
jj_disc_tpproj_..._/jj_disc_tp.{k}.h5                 # P^tp(t); per-site weighted single w*T(n,s)
jj_disc_ylmproj_..._/jj_disc_ylm.{k}.h5               # D^ylm_(ch1,ch2)(t); projected single Ttilde^(lm)(s)
```

The geometric weights are already folded in (Chunk 10), so no separate weight dump is
needed.  The notebooks only gauge-average, subtract the vacuum piece (per link/site for the
diagonal estimators from the weighted singles; per channel for `_ylmproj`), combine with the
connected output as $-C_{V,c}+C_{V,d}$, and fit/plot.

**Done in the `.ipynb` notebooks (not here):** gauge average $\langle\cdot\rangle_U$ with
jackknife; vacuum subtraction $\langle T(0)T(t)\rangle_U - \langle T(0)\rangle_U\langle T(t)\rangle_U$
(per link/site, then $w$-weighted sum for `_spproj`/`_tpproj`); assembly of the full vector
correlator $\langle JJ\rangle = -C_{V,c} + C_{V,d}$ (Eq. 3.41) by combining with the
connected output; fitting and plotting.

**Files:** the three `jj_disc_*proj_claude.cu`

---

## Interface to Existing Code

### `ConservedCurrent::apply_k` (not modified)

Template `apply_k<Gauge, LinkEl>` accepts `std::pair<int,BaseLink>` (spatial) and
`std::pair<int,Idx>` (temporal).  Canonical calling pattern from
`check_conserved_current_claude.cu`.

### `MatPoly::solve` / `from_cpu` (not modified)

```cpp
// TEMPLATE: disc_claude.cu:231-237
auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, ...);
auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch,  &D, ...);
LinOpWrapper M_DH(f_DH), M_DHD(f_DHD);
MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
MatPoly op_DHD; op_DHD.push_back(cplx(1.0), {&M_DHD});
op_DH.from_cpu<N>(DH_eta.field, eta.field);
op_DHD.solve<N>(phi.field, DH_eta.field);
```

### `GaugeExt::read` (not modified)

```cpp
// TEMPLATE: disc_claude.cu:282
U.read(dir3 + "ckpoint_lat." + std::to_string(k));
D.update(U);
```

### Geometric weights (not modified)

`kappa[il]` (`dirac_simp.h:355-360`), `link_volume[il]`, `ell[il]`, `dual_areas[s]`
(`s2n_simp.h:17,29,38`), `GetYlm(s,l,m)` after `UpdateYlm(l_max)` (`s2.h:55-56`).

---

## File Naming

All new files follow the `_claude` suffix convention:

- `src/both_3d/jj_disc_claude.cu` (single disc dump program)
- `src/both_3d/jj_conn_spproj_claude.cu` (connected, Eq. 4.29)
- `src/both_3d/jj_conn_tpproj_claude.cu` (connected, Eq. 4.32)
- `src/both_3d/jj_conn_ylmproj_claude.cu` (connected, Eq. 4.36)
- `src/both_3d/jj_disc_spproj_claude.cu` (disc post-processing, Eq. 4.29)
- `src/both_3d/jj_disc_tpproj_claude.cu` (disc post-processing, Eq. 4.32)
- `src/both_3d/jj_disc_ylmproj_claude.cu` (disc post-processing, Eq. 4.36)

---

## Open Questions

*(None blocking. Minor tuning left to implementation: number of sink times $n_s$ and number
of $\zeta$-noise hits for the diagonal projections -- a statistics/cost trade, not a design
question.)*

---

## Resolved (was Open)

- **Connected estimator (no $\gamma_5$-hermiticity).** Resolved (user spec): single source;
  source-side vector $\phi=K^{(\text{src})}D_\text{ov}^{-1}\eta$ (kernel applied after one
  base solve -> free over source insertions) and sink-side vector
  $\psi=D_\text{ov}^{-\dagger}K^{(\text{snk})\dagger}\eta$, with $C_{V,c}=\langle\psi^\dagger\phi\rangle$.
  Single bilinear -> exact, both legs $D_\text{ov}^{-1}$, no contamination.  See "Connected
  piece" above.  (Supersedes the earlier one-end / dump-$c_0$ and the sequential-source-per-
  point ideas; this organization makes the source side free and pays only sink-side solves.)
- **Connected program split.** Decided: three programs `_spproj` / `_tpproj` / `_ylmproj`
  (by projection type), since the sink-side projection folding differs.
- **Axial in scope; selected by CLI.** Decided: both axial correlators ($C_A(0\to t)$ and
  $C_A(t\to 0)$, Eqs. 3.45--3.48) are computed by the connected programs under
  `--current axial` (vector under `--current vector`).  One current per run -- **no**
  inversion reuse across currents (simplicity).  The axial formula (Q1) is confirmed correct
  and implemented **as written**.  (Supersedes the earlier "axial deferred" and "compute
  simultaneously" notes; the Sec. 4 projections are generic over current type.)
- **Diagonal projection cost (Q2).** Resolved: the `_spproj`/`_tpproj` diagonal sum uses a
  **stochastic insertion-point noise** $\zeta$ so the sink side is one composite $D^\dagger$
  solve per sink time -- avoiding the $n_\text{links}$/$n_\text{sites}$ inversion scaling.
  See "Projection dependence" and Chunk 6.
- **Disconnected vacuum subtraction.** Decided: the per-config trace dump $T^{nn'}(t)$ is
  sufficient.  A vacuum piece $\langle T\rangle_U$ may exist but is handled entirely in
  post-processing from the per-config data; no in-program subtraction.
- **Temporal coupling $\kappa^{(0)}_{t,t+1}(n)$.** Resolved: uniform in $t$, stored with a
  site index (a per-site array; candidate `lattice.dual_areas[ix]`, cf.
  `meson_pq_wall_v2_claude.cu:239`).  Exact accessor pinned when writing the
  `jj_disc_*proj_claude.cu` / `jj_conn_*proj_claude.cu` programs.  Only the projection step
  needs it.

### By Sec. 4.2/4.3

- **Per-link vs orbit-averaged output.** Resolved: the connected programs include the
  geometric weights and output the projected correlators directly; the disc dump is raw
  traces and the disc post-proc programs apply the weights (Eqs. 4.29/4.32/4.36).  In both
  cases the weighting happens in the **compiled programs**, not the notebooks.  No
  icosahedral orbit decomposition needed -- the geometric weights $A/\kappa^{(0)2}$ do the
  averaging.
- **Output size.** Resolved: raw traces are small ($O(\texttt{n\_links}\cdot Nt)$ and
  $O(\texttt{n\_sites}\cdot Nt)$ complex per config); the projected outputs are tiny.
- **Temporal currents.** Resolved: required, not deferred -- Eqs. (4.32)/(4.36) are
  entirely temporal and carry the conformal tower.  Both programs handle spatial and
  temporal links.
- **Source-timeslice averaging.** Resolved: average over $s_0$ mod $N_t$ (standard).
- **Axial correlator.** Now in scope (see decided item above); the Sec. 4 projections
  apply to all current types, so axial uses the same 4.29/4.32/4.36 machinery.
- **Production `nhits` / dilution.** Resolved: `nhits=1`, `t_block=8`, time-spin dilution,
  identical to `disc_claude.cu`.

---

## Efficiency Decisions Recorded

- **Disconnected dumps traces only.** `jj_disc_claude.cu` writes the per-config
  single-time traces $T^{nn'}(t)$, $T^{t,t+1}(n)$; the factorized correlator
  $\langle T(0)T(t)\rangle$, geometric weighting, and projections are done by the
  `jj_disc_*proj_claude.cu` post-proc programs (the notebooks only gauge-average / subtract
  vacuum / fit).  This decouples the expensive GPU solves from the cheap re-analysis and
  lets the projection scheme (Eqs. 4.29/4.32/4.36) change without re-running the GPU dump.

- **Connected: free source side, paid sink side, one current per run.** The source-side
  base solve ($D_\text{ov}^{-1}\eta$ for vector, $D_\text{ov}^{-\dagger}\eta$ for axial) is
  done once per hit; source kernels are applied after it (free over source insertion
  points/times).  The sink side costs one composite $D^\dagger$ solve per sink time
  (projection folded into the composite source: $Y_{\ell m}$ for ylm, $\zeta$-noise for the
  diagonal) -- no per-insertion inversion scaling.  The geometric weights are folded into
  the projected correlator in-program.  `--current` selects vector or axial; no inversion
  reuse across currents (simplicity).

- **Per-insertion `apply_k`.** Each `apply_k` launches `D.size-1` pole CG solves.  Batching
  over insertion points (reusing pole-$\ell$ solutions) is a future optimization; use the
  simple per-insertion approach first and profile.

- **HDF5 per-trajectory, Truncate mode**, one `.h5` per checkpoint $k$, mirroring
  `disc_claude.cu:285-286`.
