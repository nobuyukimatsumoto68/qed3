# Subtracting the contact divergence from the condensate (order parameter)

Reference: `condensate_impl_plan_claude.md` (operators), `qed3int_v2-14.pdf` Sec. 1.
Data: `data_free_vmRe*vmIm*/condensate_deter_L1/cond.h5`.

## 1. What the nonzero massless value is

The stored legs are area-weighted full traces:
$$
\texttt{etadag\_xi} = -\operatorname{tr}[A\,D_m^{-1}],\qquad
\texttt{xidag\_1mDdag\_eta} = -\overline{\operatorname{tr}[A\,(1-D_{ov})D_m^{-1}]},
$$
with $A=\operatorname{diag}(A_n)$, $V_{st}=N_t\sum_n A_n$, $\operatorname{tr}[A\mathbb 1]=2V_{st}$ (spin 2).

The condensates are
$$
\langle\sigma_{PS}\rangle = -2\,\operatorname{Re}\operatorname{tr}[A\,D_m^{-1}],\qquad
\langle\sigma_{FS}\rangle = -\operatorname{tr}[A\,D_m^{-1}] + \overline{\operatorname{tr}[A\,(1-D_{ov})D_m^{-1}]}.
$$

Free $m=0$ data: $\operatorname{Re}\operatorname{tr}[A\,D_{ov}^{-1}]/V_{st} = 1.00000000$ (8 digits). This is the
**Ginsparg-Wilson contact term**: the GW relation $\{D_{ov},\gamma\}=D_{ov}\gamma D_{ov}$ implies
$$
D_{ov}^{-1} + D_{ov}^{-\dagger} = \mathbb 1
\;\Rightarrow\;
\operatorname{Re}\operatorname{tr}[A\,D_{ov}^{-1}] = \tfrac12\operatorname{tr}[A\mathbb 1] = V_{st}\quad\text{exactly.}
$$
So $\langle\sigma_{PS}\rangle(0) = -2V_{st}$ and (using $(1-D_{ov})D_{ov}^{-1}=D_{ov}^{-1}-\mathbb 1$)
$\langle\sigma_{FS}\rangle(0) = -\operatorname{tr}[A\mathbb 1] = -2V_{st}$. Pure geometry / UV, no dynamics.
This is an ADDITIVE divergence (a constant, analytic in $m$); it is NOT removed by the $(1-D_{ov})$
flavor-symmetric combination, because FS is a different physical operator (flavor basis), not the
subtracted condensate.

## 2. The fix: GW-subtracted (half) operator

Measure the standard GW-improved bilinear (Niedermayer / Hasenfratz subtraction)
$$
\bar\psi\,\big(1-\tfrac12 D_{ov}\big)\,\psi,
\qquad
\langle\sigma^{sub}_{PS}\rangle \equiv -2\,\operatorname{Re}\operatorname{tr}\!\big[A\,(1-\tfrac12 D_{ov})D_m^{-1}\big]
= \langle\sigma_{PS}\rangle + \operatorname{Re}\operatorname{tr}[A\,D_{ov}D_m^{-1}].
$$
At $m=0$ the subtraction term is $\operatorname{tr}[A\,D_{ov}D_{ov}^{-1}]=\operatorname{tr}[A\mathbb 1]=2V_{st}$, so
$$
\langle\sigma^{sub}_{PS}\rangle(0) = -2V_{st} + 2V_{st} = 0 \quad\text{(free, massless).}
$$
The factor $\tfrac12$ is fixed by the GW contact value $\operatorname{Re}\operatorname{tr}D_{ov}^{-1}=\tfrac12\operatorname{tr}\mathbb 1$,
independent of the overlap normalization. This is the operator whose renormalization is purely
multiplicative and matches continuum $\bar\psi\psi$.

### No new solve -- build it from the legs already stored
$$
\operatorname{tr}[A\,D_{ov}D_m^{-1}] = \operatorname{tr}[A\,D_m^{-1}] - \operatorname{tr}[A\,(1-D_{ov})D_m^{-1}]
 = -\texttt{etadag\_xi} + \overline{\texttt{xidag\_1mDdag\_eta}} \cdot(\pm),
$$
i.e. the same two contractions used for PS and FS. Concretely
$$
\langle\sigma^{sub}_{PS}\rangle
= -2\,\operatorname{Re}\,\texttt{etadag\_xi}\;-\;\operatorname{Re}\!\big(\texttt{etadag\_xi}+\texttt{xidag\_1mDdag\_eta}\big)
\quad(\text{up to the stored signs}).
$$
Check on free $m=0$ data: $-2(-1608.495) + (-3216.889) = -0.10$, i.e. $-6\times10^{-5}\,V_{st}$ vs the raw
$-2V_{st}=-3217$. A factor $\sim3\times10^{4}$ smaller -- consistent with zero up to the $O(a)$ lattice
artifact at L1 (no SSB in the free theory, as required).

## 3. Complementary cross-checks / alternatives

- **Free-field subtraction.** $\langle\sigma\rangle^{ren}(m)=\langle\sigma\rangle_{\rm full}(m)-\langle\sigma\rangle_{\rm free}(m)$,
  using exactly the dense free numbers you already computed at the same $m$. Removes the additive UV piece
  at one loop; identically zero in the free limit (a built-in consistency test). Residual = genuine
  interacting condensate.
- **Mass antisymmetrization (Banks-Casher).** The divergence is even in $m$; the order parameter is the
  discontinuity:
  $$
  \Sigma = -\lim_{m\to 0^+}\lim_{V\to\infty}\tfrac12\big[\langle\sigma\rangle_{m}-\langle\sigma\rangle_{-m}\big].
  $$
  The constant + even-in-$m$ divergences cancel; needs a $-m$ run. Combine with the GW subtraction for the
  cleanest signal.

## 4. Recommendation

For the order parameter, switch the measured operator to the GW-subtracted form
$\bar\psi(1-\tfrac12 D_{ov})\psi$ (Sec. 2) -- buildable from the existing `etadag_xi` /
`xidag_1mDdag_eta` legs, zero new solves -- and verify $\to 0$ in the free massless run (done above:
$-6\times10^{-5}V_{st}$). For the parity mass $m_P$ the same construction applies with $D_m\to\tilde D_{m_P}$
and the $(1+m_P)^{-1}$ factor (the parity leg already stored as `phimm`).

## 5. References

The $(1-\tfrac12 D)$ subtracted operator (why the factor, why it is the right field):
- M. Luscher, "Exact chiral symmetry on the lattice and the Ginsparg-Wilson relation," Phys. Lett. B 428
  (1998) 342, hep-lat/9802011. The GW-compatible chiral rotation acts on $\hat\psi=(1-\tfrac12 D)\psi$, so
  the physical scalar density is $\bar\psi(1-\tfrac12 D)\psi$, not the naive $\bar\psi\psi$.
- P. Ginsparg, K. Wilson, Phys. Rev. D 25 (1982) 2649. The relation itself; source of the
  $D^{-1}+D^{-\dagger}=\mathbb 1$ contact identity used in Sec. 1.
- P. Hasenfratz, V. Laliena, F. Niedermayer, "The index theorem in QCD with a finite cutoff," Phys. Lett. B
  427 (1998) 125, hep-lat/9801021. The $\langle\bar\psi\psi\rangle$ contact term and its $\tfrac12\operatorname{tr}\mathbb 1$
  value appear explicitly.

Pedagogical write-ups (GW condensate + contact subtraction):
- C. Gattringer, C. Lang, "Quantum Chromodynamics on the Lattice," LNP 788 (Springer, 2010), Ch. 7.
- S. Chandrasekharan, U.-J. Wiese, "An introduction to chiral symmetry on the lattice," Prog. Part. Nucl.
  Phys. 53 (2004) 373, hep-lat/0405024.

Additive renormalization of $\langle\bar\psi\psi\rangle$ generally:
- M. Bochicchio, L. Maiani, G. Martinelli, G. Rossi, M. Testa, Nucl. Phys. B 262 (1985) 331. The lattice
  condensate carries a power-divergent additive piece that must be subtracted before it is an order
  parameter (Wilson-fermion context; for GW the additive divergence collapses to just the contact term).

Most directly on-point -- QED3 condensate as order parameter, with overlap, this subtraction:
- N. Karthik, R. Narayanan, "No evidence for bilinear condensate in parity-invariant three-dimensional QED
  with massless fermions," Phys. Rev. D 93 (2016) 045020, arXiv:1512.02993.
- N. Karthik, R. Narayanan, "Scale-invariance of parity-invariant three-dimensional QED," Phys. Rev. D 94
  (2016) 065026, arXiv:1606.04109 (= `obsolete/1606.04109v2.pdf` in this tree).
  Methodological note: rather than subtracting the contact term operator-side, they extract the condensate
  from the low-lying spectral density $\rho(0)$ of the massless anti-Hermitian overlap operator via a
  Banks-Casher / random-matrix analysis -- the divergence lives in the bulk/UV eigenvalues, the order
  parameter in $\rho(0)$. A third, spectral route complementary to the operator subtraction (Sec. 2) and
  the mass-antisymmetrization (Sec. 3).
