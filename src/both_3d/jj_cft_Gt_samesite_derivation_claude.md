# Analytic derivation: $G_t(t;\hat n,\hat n)$ reproduces Eq. (4.31)

Closed-form proof that the temporally projected conserved-current two-point function at coincident
spatial point, $n_1=n_2=n$, follows **exactly** from the (validated) continuum eigenbasis propagator,
reproducing qed3int_v2-11 Eq. (4.31). Equation references: `jj_cft_check_equations_claude.md`.

## 1. Inputs (all validated)

The eigenbasis propagator $G(x_1,x_2)=D_\text{cont}^{-1}$ (Eq. (C.28), code `G_prop`) satisfies, at a
coincident spatial site $\hat n$ and time separation $t>0$ (validated: same-site block is pure
$\sigma_3$ with off-diagonal machine-zero, diagonal $=$ Eq. (C.29) to $8\times10^{-16}$, and
$\hat n$-independent by homogeneity to $\sim10^{-13}$):
$$
G\big((\hat n,t),(\hat n,0)\big) = c_3(t)\,\sigma_3,
\qquad
G\big((\hat n,0),(\hat n,t)\big) = -\,c_3(t)\,\sigma_3,
$$
the sign flip on the diagonal because $G_{00}\propto\mathrm{sign}(\tau)$ is odd in $\tau=t_1-t_2$. The
coefficient is the Eq. (C.29) sum
$$
c_3(t) = \frac{1}{4\pi}\sum_{n\ge0}(n+1)\,e^{-(n+1)t}.
$$

## 2. Current and contraction

Conserved vector current, Eqs. (3.6)-(3.7):
$$
j_V^a = j_{V,+}^a + j_{V,-}^a = \eta^\dagger\sigma^a\xi - \xi^\dagger\sigma^a\eta.
$$
From the action $S_2=\eta^\dagger D\xi-\xi^\dagger D\eta$ (Eq. (1.24), $\chi$-basis (1.25)-(1.27)):
$$
\langle\xi\,\eta^\dagger\rangle = G,\qquad
\langle\eta\,\xi^\dagger\rangle = -G,\qquad
\langle\xi\,\xi^\dagger\rangle = \langle\eta\,\eta^\dagger\rangle = 0.
$$
Wick-contracting $\langle j_V^a(x_1)\,j_V^b(x_2)\rangle$ (connected): the $j_{V,+}j_{V,+}$ and
$j_{V,-}j_{V,-}$ pieces are equal; the cross pieces $j_{V,+}j_{V,-}$ vanish because they would need
$\langle\xi\xi^\dagger\rangle$ or $\langle\eta\eta^\dagger\rangle$. With the single fermion-loop minus
sign,
$$
f_{\rm cyl}^{ab}(x_1,x_2) \equiv \langle j_V^a(x_1)\,j_V^b(x_2)\rangle
 = -\,\mathcal N\,\mathrm{tr}\!\big[\sigma^a\,G(x_1,x_2)\,\sigma^b\,G(x_2,x_1)\big].
$$
The literal current gives $\mathcal N=2$ (the $++$ and $--$ terms); $\mathcal N$ is an overall
normalization absorbed into $C_J$/$C_j$ (the code `jj_cft_Gt_check_claude.cc` uses $\mathcal N=1$). The
temporal projection $e^a_3 e^b_3$ selects the $a=b=3$ component (the global time axis in this frame,
i.e. $\sigma^3$): $G_t=f^{33}$, and $G_s=f^{11}+f^{22}$.

## 3. Same-site reduction (the trace algebra)

Substitute $G(x_1,x_2)=c_3\sigma_3$, $G(x_2,x_1)=-c_3\sigma_3$:
$$
f^{ab} = -\mathcal N\,\mathrm{tr}\!\big[\sigma^a(c_3\sigma_3)\,\sigma^b(-c_3\sigma_3)\big]
       = \mathcal N\,c_3^2\,\mathrm{tr}\!\big[\sigma^a\sigma_3\,\sigma^b\sigma_3\big].
$$
Use $\sigma_3\,\sigma^b\,\sigma_3 = \eta_b\,\sigma^b$ with
$$
\eta_1=\eta_2=-1,\qquad \eta_3=+1
$$
($\sigma_3$ anticommutes with $\sigma_{1,2}$, commutes with $\sigma_3$). Then
$\mathrm{tr}[\sigma^a\sigma_3\sigma^b\sigma_3]=\eta_b\,\mathrm{tr}[\sigma^a\sigma^b]=2\,\eta_b\,\delta_{ab}$,
so
$$
f^{ab} = 2\mathcal N\,c_3^2\;\eta_a\,\delta_{ab}
       = 2\mathcal N\,c_3^2\;\mathrm{diag}(-1,-1,+1).
$$
(Check $a=b=1$: $\sigma_3\sigma_1\sigma_3=-\sigma_1$, $\mathrm{tr}[\sigma_1(-\sigma_1)]=-2=2\eta_1$. Check
$a=b=3$: $\sigma_3\sigma_3\sigma_3=\sigma_3$, $\mathrm{tr}[\sigma_3\sigma_3]=2=2\eta_3$.)

Hence
$$
G_t = f^{33} = 2\mathcal N\,c_3(t)^2,
\qquad
G_s = f^{11}+f^{22} = -4\mathcal N\,c_3(t)^2,
\qquad
\frac{G_s}{G_t} = -2 = -(D-1).
$$

## 4. The geometric sum

With $k=n+1$ and $x=e^{-t}$, $\sum_{k\ge1}k\,x^k = x/(1-x)^2$, so
$$
c_3(t) = \frac{1}{4\pi}\sum_{k\ge1}k\,e^{-kt} = \frac{1}{4\pi}\,\frac{e^{-t}}{(1-e^{-t})^2},
\qquad
c_3(t)^2 = \frac{1}{16\pi^2}\,e^{-2t}(1-e^{-t})^{-4}.
$$

## 5. Result vs Eq. (4.31)

$$
\boxed{\;
G_t(t;\hat n,\hat n) = 2\mathcal N\,c_3(t)^2
 = \frac{\mathcal N}{8\pi^2}\,e^{-2t}(1-e^{-t})^{-4}
\;}
$$
which is **exactly** Eq. (4.31),
$$
G_t(t;\hat n,\hat n) = -C_J\,e^{-\Delta t}(1-e^{-t})^{-2\Delta},
$$
with conformal dimension $\Delta=2$ and $C_J=-\mathcal N/(8\pi^2)$. Correspondingly Eq. (4.28),
$G_s=C_j(D-1)e^{-\Delta t}(1-e^{-t})^{-2\Delta}=2C_j\,e^{-2t}(1-e^{-t})^{-4}$, gives
$C_j=-\mathcal N/(8\pi^2)=C_J$, consistent with $G_s/G_t=-(D-1)=-2$.

## 6. What is exact vs convention

- **Exact (given the contraction + the $\sigma^3$ projection):** the functional form
  $G_t\propto e^{-\Delta t}(1-e^{-t})^{-2\Delta}$ with $\Delta=2$; the $\hat n$-independence (from
  $c_3$ being site-independent, validated by homogeneity); and the ratio $G_s/G_t=-2$. None of these
  depend on $\mathcal N$.
- **Convention (not predicted):** the absolute sign and magnitude of $C_J$, since
  $|C_J|=\mathcal N/(8\pi^2)$ carries the unfixed current normalization $\mathcal N$ and an overall
  sign convention (Euclidean signature / current normalization). Only the *relative* $G_s/G_t$ is
  fixed.
- **Assumptions (flagged):** (i) the contraction $f^{ab}=-\mathcal N\,\mathrm{tr}[\sigma^aG\sigma^bG]$;
  (ii) the temporal projection $e^a_3e^b_3\to\sigma^3$. Both are convention/derivation, not
  independently proven here.
- **Not covered:** the general $n_1\ne n_2$ angular structure (Eq. (4.26), the
  $\Lambda^a_\alpha(\hat n_1)\Lambda^b_\beta(\hat n_2)$ form) — that needs the off-site blocks and is
  the deferred Chunk 4; this derivation is the coincident-point ($n_1=n_2$) case only.
