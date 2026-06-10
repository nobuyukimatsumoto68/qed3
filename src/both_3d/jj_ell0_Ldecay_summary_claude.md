# The $\ell=0$ Ylm mode of $G_t$ on icosahedral lattices: why it is nonzero, and how it decays with $L$

**Summary.** The $\ell=0$ spherical-harmonic projection of the temporal current two-point
function $G_t$ vanishes exactly in the continuum (charge-conservation Ward identity). Evaluated
as a *finite-lattice quadrature over the L=1,2,4,... icosahedral points* (using the continuum,
zero-temperature propagator), it is **nonzero** — a pure angular-aliasing artifact. Its
**decay rate in cylinder time is $e^{-7\tau}$ and is independent of $L$** (fixed by icosahedral
symmetry), while its **magnitude falls as $\sim 1/L^4$** (set by the $\ell=6$ quadrature error
$E_6\sim 1/L^2$). This documents the mechanism and the numbers.

Reference: qed3int_v2-11 Sec. 4.3 (Eqs. 4.33-4.35) and App. C; validated continuum tower in
`jj_cft_ylm_check_claude.{cc,md}`; the tower implementation in `jj_local_deter_claude.cu`
(`build_Sigma_ylm`, key `ylm/l{0,1,2}`).

---

## 1. Setup

The temporal-current correlator $G_t(t;\hat n_1,\hat n_2)$ (the $\sigma_3$/"s3" channel) is
isotropic — a function of $x=\hat n_1\!\cdot\!\hat n_2$ only. Its $\ell$-projection (Eq. 4.33,
Funk-Hecke) is the Legendre coefficient

$$
\gamma_\ell(t)=\tfrac12\int_{-1}^{1} G_t(t;x)\,P_\ell(x)\,dx ,
$$

equivalently the normalized $S^2$ integral against $P_\ell(\hat n_1\!\cdot\!\hat n_2)$. The
$\ell=0$ mode is the plain sphere average,
$\gamma_0(t)=\frac{1}{4\pi}\int_{S^2}G_t\,d\Omega$.

**Continuum (zero temperature):** $\gamma_0(t)=0$ to machine precision — $G_t$ has no monopole
component. Physically this is the conserved-charge Ward identity: $\ell=0$ of the temporal
current is the total charge $Q=\int_{S^2}J^t$, which is time-independent and annihilates the
neutral vacuum, so $\langle Q(t_0)Q(t)\rangle=0$. (Confirmed in `jj_cft_ylm_check`: the 48-node
Gauss-Legendre integral gives $\gamma_0\sim10^{-16}$, $\sim10^{12}$ below $\gamma_1$.)

**Lattice estimator** (`ylm/l0`): replace $\frac{1}{4\pi}\int_{S^2}\to\sum_n A_n$ over the $L$
icosahedral dual sites, $A_n=$ `dual_areas[n]`, and use the deterministic trace
$G_t(t;\hat n_n,\hat n_{n'})=\mathrm{tr}[\sigma_3(n,t_0)\,P\,\sigma_3(n',t)\,P]$ read from the
continuum propagator $P$ (`cont_prop_L*/Dinv.0.h5`):

$$
g_0^{\rm lat}(t)=\frac{1}{(4\pi)^2}\sum_{n,n'}A_n A_{n'}\,G_t(t;\hat n_n,\hat n_{n'}).
$$

This is the same propagator as the continuum check — **only the quadrature differs** (12/42/162
lattice points vs. the exact integral).

---

## 2. Why the lattice $\ell=0$ is nonzero: aliasing

Expand $G_t(t;x)=\sum_L (2L+1)\gamma_L(t)\,P_L(x)$ and use the addition theorem
$P_L(\hat n\!\cdot\!\hat n')=\frac{4\pi}{2L+1}\sum_m Y_{Lm}(\hat n)Y_{Lm}(\hat n')$:

$$
g_0^{\rm lat}(t)=\frac{1}{4\pi}\sum_{L} \gamma_L(t)\,E_L^2,\qquad
E_L^2\equiv\sum_{m=-L}^{L}\Big|\,\sum_n A_n\,Y_{Lm}(\hat n_n)\Big|^2 .
$$

$E_L$ is the **$\ell=L$ quadrature error** of the dual-area rule (it estimates
$\int Y_{Lm}\,d\Omega=\sqrt{4\pi}\,\delta_{L0}\delta_{m0}$). So:

- $E_0^2=4\pi$ but the continuum $\gamma_0=0$ — no monopole contribution.
- If the rule integrated all $Y_{L\ge1}$ exactly, every $E_{L\ge1}=0$ and $g_0^{\rm lat}=0$.
- A finite rule has $E_{L\ge1}\neq0$ for some $L$, and those harmonics of $G_t$ **alias** into
  the $\ell=0$ bin. $g_0^{\rm lat}=\frac{1}{4\pi}\sum_{L\ge1}\gamma_L(t)E_L^2$ is exactly that
  aliasing.

---

## 3. Icosahedral symmetry: which $E_L$ survive

For any icosahedrally-symmetric set of points **with** symmetric weights, $\sum_n A_n Y_{Lm}$ is
an $I_h$-invariant vector in the spin-$L$ multiplet, hence zero **unless the $\ell=L$
representation contains the trivial rep of $I_h$** — i.e. unless $L\in\{0,6,10,12,16,\dots\}$
(the icosahedral invariant harmonics). Therefore

$$
E_L=0\ \ \text{for}\ L=1,2,3,4,5,7,8,9,11,\dots;\qquad
E_L\neq0\ \ \text{for}\ L=6,10,12,\dots
$$

This is the **spherical 5-design** property: the 12 icosahedron vertices (and every refinement,
which keeps $I_h$ symmetry) integrate all $Y_{\ell m}$ with $\ell\le5$ exactly. So the first
harmonic that leaks into $\ell=0$ is $\ell=6$, **at every refinement level $L$** — refinement
does *not* raise the design degree.

Combined with the descendant tower $\gamma_L(t)\sim e^{-(L+1)\tau}$ ($\gamma_1\sim e^{-2\tau}$,
$\gamma_2\sim e^{-3\tau}$, validated), the leading lattice $\ell=0$ is

$$
\boxed{\,g_0^{\rm lat}(t)\ \approx\ \frac{1}{4\pi}\,\gamma_6(t)\,E_6^2\ \sim\ E_6^2\,e^{-7\tau}\,}
$$

so the **decay rate is $7=6+1$, $L$-independent**, and the **magnitude is set by $E_6^2$**.

---

## 4. The geometry: $E_\ell(L)$ over the refinements

Computed directly from the lattice geometry (`pts_n{L}.dat`, `pts_dual_n{L}.dat`,
`face_dual_n{L}.dat`), with the spherical dual areas $A_n$ reproduced from
`s2n_simp.h::set_dual_areas` (L'Huilier on the dual polygons; $\sum_n A_n=4\pi$ to $10^{-12}$).
$E_\ell(L)=\big\|\sum_n A_n Y_{\ell m}(\hat n_n)\big\|$ (complex $Y_{\ell m}$):

| $L$ | $N_{\rm site}$ | $E_1\!-\!E_5$ | $E_6$ | $E_7,E_8$ | $E_6\,L^2$ | $E_6^2$ |
|----:|----:|:----:|----:|:----:|----:|----:|
| 1 | 12 | $\sim10^{-15}$ | 8.48 | $\sim10^{-15}$ | 8.48 | 71.9 |
| 2 | 42 | $\sim10^{-15}$ | 0.260 | $\sim10^{-15}$ | 1.04 | 0.0676 |
| 4 | 162 | $\sim10^{-15}$ | 0.151 | $\sim10^{-15}$ | 2.42 | 0.0228 |
| 8 | 642 | $\sim10^{-14}$ | 0.0415 | $\sim10^{-14}$ | 2.66 | 0.00172 |
| 16 | 2562 | $\sim10^{-14}$ | 0.0106 | $\sim10^{-14}$ | 2.71 | 0.000112 |
| 24 | 5762 | $\sim10^{-14}$ | 0.00473 | $\sim10^{-14}$ | 2.72 | 0.0000224 |

Readout:

- $E_{1\text{-}5}$ and $E_{7,8}$ are at machine zero for **all** $L$ — the 5-design holds at every
  refinement, and only the icosahedral invariants ($\ell=6$, then $10,12$) leak. (Confirms Sec. 3.)
- **$E_6$ decays as $\sim 2.7/L^2$:** $E_6\,L^2 \to 2.72$ for $L\ge4$. So the rule integrates the
  $\ell=6$ invariant ever better as the mesh refines (standard $O(a^2)$ quadrature convergence,
  $a\sim1/L$).
- $L=1$ is **anomalous**: $E_6=8.48$ (the 12-vertex icosahedron is a uniquely poor $\ell=6$
  integrator), then a cliff to $L=2$ and clean $1/L^2$ afterward.

---

## 5. Consequence for the $\ell=0$ correlator

Since $g_0^{\rm lat}\propto E_6^2$:

$$
\frac{g_0^{\rm lat}(L)}{g_0^{\rm lat}(L=1)} \approx \frac{E_6^2(L)}{E_6^2(1)}:\qquad
L{=}2:\ \tfrac{0.068}{71.9}\!\sim\!10^{-3},\quad
L{=}4:\ \tfrac{0.023}{71.9}\!\sim\!3\times10^{-4},\quad
L{=}8:\ \sim\!2\times10^{-5},\dots
$$

i.e. the magnitude **crashes by $\sim10^3$ from $L=1$ to $L=2$, then falls like $1/L^4$**, while
the shape stays $e^{-7\tau}$. The mode is converging to the continuum's exact zero — just
through the $\ell=6$ aliasing channel, not uniformly.

### Direct numerics at $L=1$ (continuum propagator, $t_0=0$)

$g_0^{\rm lat}(t)=\frac{1}{144}\sum_{n,n'}G_t(t;n,n')$ (equal areas $A_n=4\pi/12$ exact at $L=1$):

| $t$ ($n_t$) | $g_0^{\rm lat}$ | rate $-d\ln\lvert g_0\rvert/d\tau$ |
|----:|----:|----:|
| 1 | $-6.50\times10^{-1}$ | (transient) |
| 2 | $-3.69\times10^{-2}$ | |
| 5 | $-2.60\times10^{-4}$ | 7.41 |
| 10 | $-2.18\times10^{-7}$ | 7.02 |
| 20 | $-1.79\times10^{-13}$ | **7.00** |

vs. continuum exact-GL $\gamma_0\sim10^{-16}$ (machine zero). The rate cleanly $\to 7$, as predicted.

---

## 6. Practical notes

- **Rate vs. magnitude.** The $e^{-7\tau}$ rate is locked by icosahedral symmetry and will *not*
  move with $L$ (or with more statistics) — it is the signature of the $\ell=6$ leak. Only the
  magnitude converges ($\sim1/L^4$).
- **Watch the quadrature weights.** The decay only appears with the *real* dual areas $A_n$
  (`dual_areas[n]`). A naive equal-area shortcut ($A_n=4\pi/12$) is exact only at $L=1$; reusing
  it (or a stale $1/(4\pi)^2$ normalization / fixed site count) at $L\ge2$ will mask the
  $E_6^2\!\sim\!1/L^4$ fall and can make the magnitude look $L$-independent. Use the geometry's
  `dual_areas` and the correct $\frac{1}{(4\pi)^2}\sum_{n,n'}A_nA_{n'}$ prefactor.
- **Larger $L$.** Helps (it is already shrinking by $L=2$); the bottleneck is the $\ell=6$
  invariant, whose error falls as $1/L^2$. There is no design-degree improvement to wait for —
  the gain is purely the smooth $O(a^2)$ convergence of the $\ell=6$ integration.

---

## 7. Reproducibility

- Geometry / $E_\ell(L)$: load `pts_n{L}.dat` (sites), `pts_dual_n{L}.dat` (dual sites),
  `face_dual_n{L}.dat` (dual polygons) from `../../geometry/data/`; build $A_n$ via L'Huilier per
  `set_dual_areas`; $E_\ell=\sqrt{\sum_m|\sum_n A_n Y_{\ell m}(\hat n_n)|^2}$ with
  `scipy.special.sph_harm_y(l,m,polar,azimuth)`. (One Python script, $L=1\!-\!24$, seconds.)
- $g_0^{\rm lat}(t)$ at $L=1$: from `cont_prop_L1/Dinv.0.h5` (`Dm_inv/{real,imag}`, row-major
  $N\times N$, $N=N_t\,n_{\rm site}\,N_s=128\cdot12\cdot2$), $t_0=0$, the $\sigma_3$ trace above.
- The compiled tower is `jj_local_deter_claude.cu` key `h0/t0_b/ylm/l0/Vpp` (uses real
  `dual_areas`, double sum) — should reproduce the $E_6^2$ scaling across $L$.

---

## 8. Observed in `jj_local_deter` (continuum prop) — CONFIRMED, and the earlier "not observed" cause

The compiled tower `corr_deter_local_cont_L{1,2,4}/corr.0.h5` (key `ylm/l0/Vpp`, real `dual_areas`,
diagonal-$m$ $g_\ell$) **does** show the decay. Measured $|g_0|$ ($t_0=0$):

| $t$ ($n_t$) | $L=1$ | $L=2$ | $L=4$ |
|----:|----:|----:|----:|
| 5  | $2.60\times10^{-4}$ | $2.79\times10^{-6}$ | $8.13\times10^{-8}$ |
| 10 | $2.18\times10^{-7}$ | $2.24\times10^{-10}$ | $6.92\times10^{-11}$ |

and $g_0/g_1$ at $dt=5$: $0.282 \to 4.6\times10^{-3}\to1.3\times10^{-4}$ — the magnitude crashes
$\sim10^2$ from $L=1\to2$ and keeps falling, exactly the $E_6^2$ aliasing channel (Sec. 5). The
$\ell=1,2$ towers stay put ($|g_1|,|g_2|$ at $dt=5$ are $\sim9\times10^{-4}$ for all $L$). Confirms
the prediction; the data was always correct.

**Why it looked $L$-independent ("we did not observe this").** A *notebook plotting bug*, not the
data: the Ylm-tower overlay (`checkCFT_claude.ipynb`, cell `99d8cb78`) reloaded `g1,g2` per $L$ inside
the loop but plotted `np.abs(g0[dt])` using the **global** `g0` left over from the $L=1$ single-$L$
cell — so $\ell=0$ was drawn as the $L=1$ value for every $L$. Fixed by loading `g0` per $L$, plus a
dedicated $\ell=0$-decay cell (`cft_ylm_l0_Ldecay_cont.pdf`: $|g_0|$ vs $t$ for $L=1,2,4$ with the
$e^{-7t}$ reference, and the magnitude-ratio printout). No recomputation needed — the `--out-tag cont`
h5 files already hold the result, separate from the lattice-prop `corr_deter_local_L*`.
