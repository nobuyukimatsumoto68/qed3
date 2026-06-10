# Plan + equations — check Appendix C identities and the Eq. (4.35) descendant tower

Source: qed3int_v2-11 Appendix C (real spherical harmonics) and Sec. 4.3 (Eqs. 4.33-4.35).
Builds on the validated $G_t$ machinery (`jj_cft_fcyl_check_claude.cc`, `jj_cft_sptp_check_claude.cc`):
$G_t(t;\hat n_1,\hat n_2)=f^{33}=-\mathrm{tr}[\sigma^3 G(x_1,x_2)\sigma^3 G(x_2,x_1)]$, verified isotropic
(depends only on $\hat n_{12}$) and equal to $\hat n_1\!\cdot T\hat n_2$ to machine precision.

## Appendix C equations (to check)

Real spherical harmonics without Condon-Shortley phase, Eq. (C.1):
$$
Y_{\ell,m}(\theta,\phi)=\sqrt{\tfrac{2\ell+1}{4\pi}\tfrac{(\ell-|m|)!}{(\ell+|m|)!}}\,P_\ell^{|m|}(\cos\theta)\times
\begin{cases}\sqrt2\cos|m|\phi & m>0\\ 1 & m=0\\ \sqrt2\sin|m|\phi & m<0\end{cases}
$$
with $P_\ell^{|m|}$ the associated Legendre functions (C.2)-(C.3). Identities to verify:
- (C.5) orthonormality: $\int d^2 n\,Y_{\ell_1 m_1}Y_{\ell_2 m_2}=\delta_{\ell_1\ell_2}\delta_{m_1 m_2}$.
- (C.6) $\int d^2 n\,Y_{\ell m}=\sqrt{4\pi}\,\delta_{\ell0}\delta_{m0}$.
- moment identities (C.7)/(C.9): $\int d^2n\,Y_{\ell m}n^a$ (spin-1, $\propto\delta_{\ell1}$) and
  $\int d^2n\,Y_{\ell m}n^a n^b$ (spin-0+spin-2, $\propto\delta_{\ell0}\delta^{ab}+\delta_{\ell2}\dots$).
- **selection rules (C.15)-(C.17)** (the ones feeding 4.35):
$$
\tfrac1{4\pi}\!\int\! Y_{\ell_1 m_1}(n_1)Y_{\ell_2 m_2}(n_2)\,\{1,\ \hat n_{12},\ \hat n_{12}^2\}
=\Big\{\delta_{\ell_1 0}\delta_{\ell_2 0}\delta_{m_1 0}\delta_{m_2 0},\ \ \tfrac13\delta_{\ell_1 1}\delta_{\ell_2 1}\delta_{m_1 m_2},\ \
\tfrac13\delta_{\ell_1 0}\delta_{\ell_2 0}+\tfrac{2}{15}\delta_{\ell_1 2}\delta_{\ell_2 2}\delta_{m_1 m_2}\Big\}.
$$
These follow from the factorized single-sphere moments: $\int\!\!\int YY\hat n_{12}^k=\sum_{a\dots}(\int Y n^{a\dots})(\int Y n^{a\dots})$,
so we verify them via the Cartesian moments $\int Y$, $\int Y n^a$, $\int Y n^a n^b$ (no need for the
frame-specific $I_1,I_2$ maps, since $\hat n_{12}=\sum_a n_1^a n_2^a$ is basis-independent).

## Eq. (4.34)/(4.35) tower

$\Delta=2$, projecting $G_t$ onto $Y_{\ell m}$ (4.33). By Funk-Hecke + isotropy of $G_t$:
$$
G^t_{\ell_1 m_1;\ell_2 m_2}(t)=\delta_{\ell_1\ell_2}\delta_{m_1 m_2}\,G^t_{\ell\ell}(t),\qquad
G^t_{\ell\ell}(t)=\tfrac12\!\int_{-1}^1 G_t(t;x)P_\ell(x)\,dx .
$$
Prediction (4.35), with our $C_j=-1/(8\pi^2)$:
$$
G^t_{11}(t)\to -\tfrac{C_j}{3}e^{-2t}=\tfrac{1}{24\pi^2}e^{-2t},\qquad
G^t_{22}(t)\to -\tfrac{4C_j}{5}e^{-3t}=\tfrac{1}{10\pi^2}e^{-3t},\qquad
G^t_{00}\to 0\ (\text{coeff }2(2-\Delta)/3=0),
$$
with $\ell=3$ at $O(e^{-4t})$. Convention-free checks: decay rates $\to(2,3)$ for $\ell=(1,2)$;
ratio $[G^t_{22}e^{3t}]/[G^t_{11}e^{2t}]\to 12/5=2.4$; $G^t_{00}/G^t_{11}\to0$.

## Files / chunks

- NEW `jj_cft_ylm_check_claude.cc`:
  - Part A: real $Y_{\ell m}$ (GSL `sphPlm`); verify (C.5) and (C.15)-(C.17) by sphere quadrature
    (Gauss-Legendre in $\cos\theta$ x uniform $\phi$) via the single-sphere Cartesian moments.
  - Part B: $g(t,x)=G_t$ from the propagator at Gauss-Legendre nodes $x_k$ (equatorial pairs,
    overflow-safe); $G^t_{\ell\ell}(t)=\tfrac12\sum_k w_k g(t,x_k)P_\ell(x_k)$ for $\ell=0..3$ at IR
    $t=2,3,4,5$; report values, effective decay rates, $G^t_{11}e^{2t}$, $G^t_{22}e^{3t}$, ratio.
- NEW `run_jj_cft_ylm_check_claude.sh` (`g++ -x c++` + GSL + Boost; NO CUDA).

## Numerics / assumptions

- nmax=40, equatorial pairs (z=0) for the propagator (IR / overflow-safe). GL with N=48 in $x$;
  sphere quadrature 32 ($\cos\theta$) x 32 ($\phi$). $\Delta=2$, $C_j=-1/(8\pi^2)$ (our $\mathcal N=1$).
- The $\delta_{\ell_1\ell_2}\delta_{m_1 m_2}$ structure of the tower is guaranteed by Funk-Hecke +
  the verified isotropy; Part A independently confirms the Appendix C selection rules that the
  paper's derivation uses.
