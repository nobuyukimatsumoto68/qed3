# JJ current-correlator CFT check — relevant equations (qed3int_v2-11.pdf)

Records the equations needed to check the Section 4 CFT formulas for the conserved-current
two-point function using the (now validated) continuum eigenbasis propagator. First target:
**Eq. (4.31)** for $G_t$ (defined in Eq. (4.30)).

## Action and field content (Sec. 1.2.1)

Continuum action, Eq. (1.24):
$$S_2 = \eta^\dagger D_\text{cont}\,\xi - \xi^\dagger D_\text{cont}\,\eta.$$

Flavor-symmetric variables, Eqs. (1.25)-(1.27):
$$\chi_1\equiv\xi,\quad \chi_2\equiv\eta,\qquad \bar\chi_1\equiv\eta^\dagger,\quad \bar\chi_2\equiv-\xi^\dagger,$$
$$S_2 = \bar\chi_1 D_\text{cont}\chi_1 + \bar\chi_2 D_\text{cont}\chi_2.$$

**Propagators.** With $\langle\chi_i\bar\chi_j\rangle=\delta_{ij}D_\text{cont}^{-1}\equiv\delta_{ij}G$:
$$\langle\xi\,\eta^\dagger\rangle = G,\qquad \langle\eta\,\xi^\dagger\rangle = -G,\qquad \langle\xi\,\xi^\dagger\rangle=\langle\eta\,\eta^\dagger\rangle=0.$$
Here $G=G(x_1,x_2)$ is exactly the eigenbasis propagator computed by `cont_prop_eigbasis_claude.cpp`
(`G_prop`), $G=D_\text{cont}^{-1}$. The $0$'s are because the action is flavor-diagonal in $\chi$.

## Conserved currents (Sec. 3.1, Eqs. (3.6)-(3.11))

$\sigma^a$ ($a=1,2,3$) are the Pauli matrices (the 2-component $\gamma^a$). $T_r$ are flavor generators.

1. Vector $U(1)^2$, Eqs. (3.6)-(3.7):
$$j_{V,\pm}^a \equiv \bar\Psi_f\gamma^a\frac{1\pm\gamma_{4,5}}{2}\Psi_f
 = \begin{cases} \eta_f^\dagger\sigma^a\xi_f & (+)\\ -\xi_f^\dagger\sigma^a\eta_f & (-)\end{cases}.$$

2. Axial $\Gamma=\gamma_4,\gamma_5$, Eqs. (3.8)-(3.9):
$$j_{A,\Gamma}^a \equiv \bar\Psi_f\gamma^a\Gamma\Psi_f
 = \begin{cases} \eta_f^\dagger\sigma^a\eta_f-\xi_f^\dagger\sigma^a\xi_f & (\Gamma=\gamma_4)\\
                 -i(\eta_f^\dagger\sigma^a\eta_f+\xi_f^\dagger\sigma^a\xi_f) & (\Gamma=\gamma_5)\end{cases},$$
with $j_{A,\pm}^a\equiv\tfrac12(j_{A,\gamma_4}^a\pm i j_{A,\gamma_5}^a)$.

3. Flavor $SU(N_f/2)^2$, Eqs. (3.10)-(3.11):
$$j_{r,\pm}^a \equiv \bar\Psi_f\gamma^a\frac{1\pm\gamma_{4,5}}{2}(T_r)_{fg}\Psi_g
 = \begin{cases} \eta_f^\dagger\sigma^a(T_r)_{fg}\xi_g & (+)\\ -\xi_f^\dagger\sigma^a(T_r)_{fg}\eta_g & (-)\end{cases}.$$

**Current used for $G_t/G_s$:** the conserved **vector** current
$j_V^a = j_{V,+}^a + j_{V,-}^a = \eta^\dagger\sigma^a\xi - \xi^\dagger\sigma^a\eta$ (single flavor pair; the
gauge sector couples to this $U(1)$).

## Cylinder current correlator (Eq. (4.26))

With $t\equiv t_1-t_2$, $\hat n_{12}\equiv\hat n_1\cdot\hat n_2$, $\Lambda^a_\alpha(\hat n)$ the local
frame (Eq. (4.25), rows $\hat\theta,\hat\phi,\hat n$), the conserved-current two-point function on the
cylinder is, Eq. (4.26):
$$f_{\rm cyl}^{ab}(t;\hat n_1,\hat n_2) = C_j\,\Lambda^a_\alpha(\hat n_1)\Lambda^b_\beta(\hat n_2)\,
 e^{-\Delta t}\Big[\delta^{\alpha\beta} - \frac{2(\hat n_1-e^{-t}\hat n_2)^\alpha(\hat n_1-e^{-t}\hat n_2)^\beta}
 {1-2\hat n_{12}e^{-t}+e^{-2t}}\Big](1-2\hat n_{12}e^{-t}+e^{-2t})^{-\Delta}.$$
The conserved current has dimension $\Delta = D-1 = 2$ in $D=3$.

## Spatially and temporally projected correlators (Sec. 4.2-4.3)

Spatial projection, Eqs. (4.27)-(4.28):
$$G_s(t;\hat n_1,\hat n_2)\equiv(\delta^{ab}-e^a_3 e^b_3)f_{\rm cyl}^{ab},\qquad
  G_s(t;\hat n,\hat n) = C_j (D-1)\,e^{-\Delta t}(1-e^{-t})^{-2\Delta}.$$

**Temporal projection (the target), Eqs. (4.30)-(4.31):**
$$G_t(t;\hat n_1,\hat n_2)\equiv e^a_3 e^b_3\,f_{\rm cyl}^{ab},\qquad
  \boxed{\,G_t(t;\hat n,\hat n) = -C_J\,e^{-\Delta t}(1-e^{-t})^{-2\Delta}\,}.$$
Here $e^a_3$ is the local-frame unit vector along the radial / cylinder-time direction (the $\hat n$
row of $\Lambda$, Eq. (4.22)); $e^a_3 e^b_3 f^{ab}$ selects the $\hat n$-$\hat n$ ("33") component.

Real-spherical-harmonic descendants (for later, Eqs. (4.33)-(4.35)) confirm $\Delta=2$:
$G^t_{\ell\ell}\sim C_j[-\tfrac13\delta_{\ell1}e^{-2t} - \tfrac45\delta_{\ell2}e^{-3t}+O(e^{-4t})]$.

## Contraction in terms of the propagator $G$

Wick-contracting $\langle j_V^a(x_1)\,j_V^b(x_2)\rangle$ with the propagators above (the $++$ and $--$
pieces are equal, the $+-$ cross terms vanish since $\langle\xi\xi^\dagger\rangle=\langle\eta\eta^\dagger\rangle=0$):
$$f_{\rm cyl}^{ab}(x_1,x_2) = -\,\mathcal N\,\mathrm{tr}\big[\sigma^a\,G(x_1,x_2)\,\sigma^b\,G(x_2,x_1)\big],
  \qquad \mathcal N>0\ \text{absorbed into}\ C_J.$$
(Fermion-loop minus sign; $\mathcal N=2$ from $++$ plus $--$. The overall constant, including its
sign, is not predicted — it is absorbed into $C_J$/$C_j$; the physics checks are the $t$-shape, the
$\hat n$-independence, and the ratio $G_s/G_t=-(D-1)=-2$.)

## Analytic reduction at $n_1=n_2=n$ (sanity target)

The same-site time-separated block is pure $\sigma_3$ (validated): $G((\hat n,t),(\hat n,0))=c_3(t)\sigma_3$,
$G((\hat n,0),(\hat n,t))=c_3(-t)\sigma_3=-c_3(t)\sigma_3$, with (from Eq. (C.29))
$$c_3(t)=\frac{1}{4\pi}\sum_{n\ge0}(n+1)e^{-(n+1)t}=\frac{1}{4\pi}\frac{e^{-t}}{(1-e^{-t})^2}.$$
Then, with $e^a_3 e^b_3\to\sigma^3$ (the global time axis in our frame):
$$G_t = f_{\rm cyl}^{33} = -\mathrm{tr}[\sigma^3(c_3\sigma_3)\sigma^3(-c_3\sigma_3)] = 2c_3(t)^2
 = \frac{1}{8\pi^2}\,e^{-2t}(1-e^{-t})^{-4},$$
matching Eq. (4.31) with $\Delta=2$, $C_J=-1/(8\pi^2)$. Likewise
$f^{11}=f^{22}=-\mathrm{tr}[\sigma^1(c_3\sigma_3)\sigma^1(-c_3\sigma_3)]=-2c_3^2$, so
$G_s=f^{11}+f^{22}=-4c_3^2$ and $G_s/G_t=-2$. The `.cc` computes $f^{ab}$ numerically from `G_prop`
(no use of the pure-$\sigma_3$ shortcut) and compares to these.

## Conventions / identifications used in the `.cc`

- $t$ = cylinder (radial-quantization) time = the $\tau$ of `G_prop` (dimensionless); $\Delta=2$.
- Propagator: `G_prop` (eigenbasis, Boost $\xi$), normalization $1/(4\pi c^2_{|m|,n})$ — reproduces
  $c_3(t)$ above and Eq. (C.29) exactly.
- Temporal projection $e^a_3 e^b_3$ realized as the $a=b=3$ ($\sigma^3$) component, justified by the
  same-site block being pure $\sigma_3$ in this frame (global time axis). [CONFIRM]
- Overall constant $C_J$ (including sign) is convention; checks are shape, $\hat n$-independence,
  $G_s/G_t=-2$. [CONFIRM]
