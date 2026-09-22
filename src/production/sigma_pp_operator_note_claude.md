# Back-to-back two-meson operator $\sigma(p)\sigma(-p)$ on $S^2\times\mathbb{R}$

Goal: a two-meson $0^{++}$ operator that couples to the genuine two-particle sector with **no**
vacuum overlap and **suppressed** ground-single-meson overlap, complementing the zero-momentum
$\sigma_{00}^2$ that is contaminated by the single-meson triangle.

## 1. Why $\sigma_{00}^2$ alone fails (recap)

Both $\sigma_{00}$ and $\sigma_{00}^2$ are built from the same $l=0$ (rotationally invariant)
vertex. $\sigma_{00}^2$ is a genuine 4-fermion (two-meson) operator, but:

- it has a nonzero vacuum piece $\langle\sigma_{00}^2\rangle = 2(D_S^2+D_S')\neq 0$, and
- it overlaps the ground single meson through the cubic triangle
  $\langle\sigma_{00}\,\sigma_{00}\,\sigma_{00}\rangle\neq 0$ (diagram A / $S_S$).

On $S^2$ "momentum" is angular momentum $Y_{lm}$. The $l=0$ single-meson tower
($2E_0=m_\sigma,\;E_0+E_1,\dots$) is dense and reaches up through $2m_\sigma$, so the two-meson
level stays buried no matter how much we Lanczos-inflate the $l=0$ operators.

## 2. The operator

Angular-momentum-projected scalar density (contact-subtracted $\sigma=\bar\psi S\psi-\tfrac12$):

$$
\sigma_{lm}(t) = \sum_x w_x\, Y_{lm}^*(x)\,\sigma(x,t),
\qquad w_x = \text{dual area}(x).
$$

Couple two of them to total spatial $L=0$ (the $0^{++}$ singlet of two angular momenta $l$):

$$
O^{(2)}_l(t) = \frac{1}{\sqrt{2l+1}}\sum_{m=-l}^{l} (-1)^{l-m}\,
                 \sigma_{l,m}(t)\,\sigma_{l,-m}(t).
$$

Quantum numbers: parity $(-1)^l(-1)^l=+1$, total $L=0$, Bose-symmetric (coupling two identical
integer $l$ to $L=0$ has exchange sign $(-1)^{2l}=+1$). So $O^{(2)}_l$ stays in the $0^{++}$ sector
for every $l$. The $l=0$ member is exactly $O^{(2)}_0=\sigma_{00}^2$.

## 3. What $l>0$ actually buys (honest accounting)

$O^{(2)}_l$ is **not** free of the vacuum or the single meson -- both leak in:

- **It has a vacuum constant.** $\langle O^{(2)}_l\rangle=\sum_m(-1)^{l-m}
  \langle\sigma_{l,m}\sigma_{l,-m}\rangle$. The *disconnected* factor $\langle\sigma_{l,m}\rangle=0$
  (odd projector, $\sum_x w_x Y_{lm}=0$), but the **connected** equal-time piece -- the
  $l$-component of the meson two-point -- is nonzero. So like $\sigma_{00}^2$ it needs the identity
  operator / vacuum handling.
- **It still overlaps the ground single meson.** The triangle $\langle\sigma_{00}\,O^{(2)}_l\rangle
  \propto\sum_m(-1)^{l-m}\langle\sigma_{00}\,\sigma_{l,m}\,\sigma_{l,-m}\rangle$. For $l=1$ the
  angular integral $\int Y_{00}\,n_\mu n_\nu\propto\delta_{\mu\nu}$ is nonzero: the trace part
  carries $l=0$ and reaches the $s$-wave ground. (Typically *reduced* vs
  $\langle\sigma_{00}^3\rangle$, not absent.)

The genuine value is variational:

- **A linearly-independent $0^{++}$ two-meson interpolator** with *different* overlaps -- strong on
  two-meson scattering states of nonzero relative angular momentum, differently weighted on the
  single-meson tower than $\sigma_{00}^2$. In a GEVP the extra independent direction lets the
  generalized eigenvectors separate the near-degenerate two-meson level from the dense single tower
  that buried it in the $\{1,\sigma_{00},\sigma_{00}^2\}+$Lanczos runs. No single operator "avoids"
  a state; the basis diagonalizes them all.

## 4. Real $l=1$ Cartesian form + a general mode-space Wick engine

For $l=1$ the $L=0$ singlet of two real vectors is just the **dot product** over the three
Cartesian components $\mu\in\{x,y,z\}$ (no complex $Y_{1m}$, no CG signs -- manifestly real and
rotation invariant):

$$
\sigma_{1,\mu}(t) = \sum_x w_x\, n_\mu(x)\,\sigma(x,t),
\qquad
O^{(2)}_1(t) = \sum_{\mu=x,y,z}\sigma_{1,\mu}(t)\,\sigma_{1,\mu}(t),
$$

with $n_\mu(x)$ the Cartesian unit-sphere coordinate (site position, $|n|=1$), and vertex

$$
W_{1\mu} = \operatorname{diag}_x\!\big(w_x\, n_\mu(x)\big)\otimes\mathbf{1}_\text{spin},
\qquad
\Phi_{1\mu}(t) = V(t)^\dagger\,W_{1\mu}\,V(t).
$$

Because `diags_pair` bakes in one vertex per timeslice, the two distinct legs $\Phi_{1x},\Phi_{1y}$
do not fit it. Instead use the **general mode-space Wick sum** (the mode-space analogue of the
validated `wick4` in `point_source_full_claude.py`): for $n$ bilinears with vertices $\Phi_i$ at
times $t_i$,

$$
\Big\langle\textstyle\prod_i \bar\psi\Phi_i\psi\Big\rangle
  = \sum_{\pi\in S_n}\ \prod_{\text{cycles }c}(-1)\,
    \operatorname{Tr}\!\Big[\textstyle\prod_{i\in c}\Phi_i\,\tilde\tau(t_i,t_{\pi(i)})\Big],
\qquad
\tilde\tau(a,b)=\tau(a,b)-\tfrac12\,\delta_{ab}\,\mathbf 1 .
$$

$\tilde\tau$ carries the contact subtraction on every equal-time pair (self- and cross-loop), exactly
as $\tilde\tau$ does now. Setting $G_{ij}=\Phi_i\,\tilde\tau(t_i,t_j)$ this is `wick4`/`wick3`
verbatim. **Cross-check:** `wick4` on $(\Phi_{00}(t),\Phi_{00}(t);\Phi_{00}(0),\Phi_{00}(0))$ must
reproduce $2\sum_i W^{10}_i\,\text{diags\_pair}_i$.

GEVP blocks needed (source time $s$, sink $t=s+dt$, translation-averaged):

- $C_{\sigma\sigma}$ = wick2 (single meson) -- already have it as `C11`.
- $C_{\sigma,\sigma^2}$ = wick3 triangle $2C'$ -- already have it as `C12`.
- $C_{\sigma^2\sigma^2}$ = wick4 -- already have it as `C22`.
- $C_{\sigma,O_1}=\sum_\mu$ wick3$\big(\Phi_{00}(t);\Phi_{1\mu}(s),\Phi_{1\mu}(s)\big)$ -- NEW.
- $C_{\sigma^2,O_1}=\sum_\mu$ wick4$\big(\Phi_{00}(t)^2;\Phi_{1\mu}(s)^2\big)$ -- NEW.
- $C_{O_1O_1}=\sum_{\mu\nu}$ wick4$\big(\Phi_{1\mu}(t)^2;\Phi_{1\nu}(s)^2\big)$ -- NEW.
- identity row: $\langle\sigma_{00}\rangle=0$, $\langle\sigma_{00}^2\rangle$ and $\langle O_1\rangle$
  measured as equal-time one-points $\langle\sigma_a\sigma_b\rangle=-\operatorname{Tr}[\Phi_a\tilde\tau\,\Phi_b\tilde\tau]$
  (the $\operatorname{Tr}[\Phi\tilde\tau]$ tadpole vanishes).

## 5. Implementation + expectation

Basis $\{\,1,\ \sigma_{00},\ \sigma_{00}^2,\ O^{(2)}_1\,\}$ ($4\times4$), identity carrying the vacuum
(its off-diagonals are the measured one-points above). Script `two_meson_gevp_pp_free_claude.py`:

1. `nmu = load_vec3(pts_n{L})` Cartesian coords -> `W1mu = repeat(w_x n_mu, NS)` (3 vertices).
2. general `wickn(Phis, times, tau)` in mode space (n=2,3,4).
3. assemble the 6 blocks + identity, solve the $4\times4$ GEVP vs $t_0$.

Free-field expectation: level 0 = vacuum ($=0$), level 1 $=m_\sigma$, and -- the test -- whether a
two-meson level near $2m_\sigma$ ($=0.786$ at L2) now separates from the excited-single band. This is
the diagnostic for whether the relative-$l=1$ interpolator resolves what Lanczos alone could not; if
$l=1$ is insufficient, add $l=2$ ($O^{(2)}_2$, same engine).

## Reference
- variational two-particle operators / single-particle contamination: Luscher-Wolff (1990);
  Blossier et al. arXiv:0902.1265; Hadron Spectrum Collab.
- distillation vertex projection: Peardon et al. 0905.2160.
