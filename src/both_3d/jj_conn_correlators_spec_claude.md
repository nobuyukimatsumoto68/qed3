# `jj_conn_correlators_claude.cu` -- clean specification (equation <-> code)

Clean, self-contained spec of the **unified connected** current-current correlator program, written so each
estimator equation can be checked **line-to-line** against the code.  Reference: `qed3int_v2-10.pdf`
(equation numbers below) and `conserved_current_correlators_impl_plan_v3_claude.md` (design history).

Disconnected is a **separate** program (`jj_disc_claude.cu`) -- not here.

Conventions: two-component spinors ($\xi/\eta$); overlap $D_\text{ov}$ (Zolotarev); $Z_2$ stochastic source
$\eta$ with $E[\eta\eta^\dagger]=\mathbb{1}$; $N=N_s\,N_\text{site}\,N_t$.  All overlap applies/solves go through
the multishift `_ms` entry points; the kernel $K$ is the multishift `apply_k_ms` (via `op_K.from_cpu`).

---

## 1. Operators  (code `:240-284`)

| symbol | definition | code |
|---|---|---|
| $D_\text{ov}$ | massless overlap | `Fermion D(DW, 0, 21)` |
| $D_m$ | $D_\text{ov}+m$ (Eq. 3.60) | `Fermion Dm(DW, valence_mass, 21)` |
| $\tilde D_{m_P}$ | $D_\text{ov}+m/(1-m)$ (Eq. 3.63) | `Fermion Dtil(DW, m/(1-m), 21)` |
| $1-D_\text{ov}$ | axial GW factor | `op_oneMinusD` (`:265`) |
| $K$ | conserved-current kernel | `kop` + `op_K`; `kop.set_temporal(U,t,n,dag)` / `set_spatial(U,t,lk,dag)` |

**Solve conventions** (uniform op-set, `:248-284`).  For overlap $X\in\{D_\text{ov},D_m,\tilde D_{m_P}\}$ the
pieces are `op_X` $=X$, `op_XH` $=X^\dagger$, `op_Xsq` $=X^\dagger X$.  Because $X^\dagger X$ is solved by CG:
$$
X^{-1}b=\texttt{op\_Xsq.solve}(\texttt{op\_XH}(b))=(X^\dagger X)^{-1}X^\dagger b,\qquad
X^{-\dagger}b=\texttt{op\_Xsq.solve}(\texttt{op\_X}(b))=(X X^\dagger)^{-1}X\,b .
$$
So a **forward** solve $X^{-1}$ uses the `op_XH` RHS-former; a **dagger** solve $X^{-\dagger}$ uses `op_X`.

---

## 2. Weights  (code `:278-301`)

$$
w^{tp}_n=\frac{A_n}{\kappa_t(n)^2}\ \text{(Eq. 4.32)},\qquad
w^{sp}_{il}=\frac{A_{nn'}}{\kappa^{(0)}_{nn'}{}^{2}}\ \text{(Eq. 4.29)},\qquad
W_\ell(n)=\frac{A_n}{\kappa_t(n)}\sum_{m=-\ell}^{\ell}Y_{\ell m}(\hat n)\ \text{(Sec. 3.6, }\kappa^1).
$$
Code: `w_tp[n]=dual_areas[n]/kappa_t[n]^2`, `w_sp[il]=link_volume[il]/kappa[il]^2`,
`W_ell[l][n]=dual_areas[n]*(sum_m Ylm_real)/kappa_t[n]`.  $\ell=0,1,2$ (`N_ELL=3`); $\ell=0\approx0$ (charge
conservation, a check).  $1/(4\pi)$ folded at output (`inv4pi`, in `write_corr`).

---

## 3. Shared forward leg  (code `:381-384`)

Per hit: $\eta=$ `fill_z2_source`; then the leg shared by **all** vector projections,
$$
\phi'=D_m^{-1}\eta:\quad \texttt{op\_DmH.from\_cpu(tmp,eta)};\ \texttt{op\_Dmsq.solve(phi,tmp)} .
$$
`flavor` $=$ ($m$ real $\neq0$); `parity` $=$ ($m$ imaginary $\neq0$).

---

## 4. Vector estimators

### 4.1 Temporal tp  (++ `:387-428`, parity (--) `:430-465`)

Insertion-diagonal over sites $n$.  The connected trace (Eqs. 3.52, 3.64):
$$
G^{t,++}(t_0\!\to\! t)=\frac1{4\pi}\sum_n w^{tp}_n\,\psi^{tp}_{n}(t_0)^{\dagger}\,K(n,t)\,\phi',\qquad
\psi^{tp}_n(t_0)=D_m^{-\dagger}K^\dagger(n,t_0)\eta .
$$
Code: source leg (`:392-396`) `set_temporal(U,t0,n,dag=true)` -> `op_K(rho,eta)` $[=K^\dagger(n,t_0)\eta]$ ->
`op_Dm(tmp,rho)` -> `op_Dmsq.solve(psi_tp[ITP(n,b)],tmp)` $[=D_m^{-\dagger}\cdot]$.  Sink (`:405-419`)
`set_temporal(U,t,n,dag=false)` -> `op_K(kphi,phi)` $[=K(n,t)\phi']$ -> `Ctp[b][(t-t0_b)\bmod N_t] += w_tp[n]*psi_tp[ITP(n,b)].dag(kphi)`.

$(--)$: massless/$m_F$ -> $G^{--}=\overline{G^{++}}$ (`write_corr_conj`).  Parity -> independent tilde
(`:430-465`): $\phi'_{-}=\tilde D_{m_P}^{-\dagger}\eta$ (`phimm`), $\psi=\tilde D_{m_P}^{-1}K(n,t_0)\eta$ (forward),
sink $K^\dagger(n,t)\,\phi'_{-}$.

**Sign:** $\mathrm{Re}\,G_t<0$ (Eq. 4.31).

### 4.2 Spherical-harmonic tower ylm  (++ within `:387-428`, parity within `:430-465`)

Keep $\ell_1=\ell_2=\ell$, sum $m_1,m_2$ -> tower; the $m$-sum collapses **before** the solve (Sec. 3.6):
$$
G^{++}_\ell(t)=\frac1{4\pi}\,\Psi_\ell(t_0)^\dagger\,\Phi_\ell(t),\quad
\Psi_\ell(t_0)=D_m^{-\dagger}\!\sum_n W_\ell(n)K^\dagger(n,t_0)\eta,\quad
\Phi_\ell(t)=\sum_n W_\ell(n)K(n,t)\phi' .
$$
**Shared with tp:** the SAME `rho`$=K^\dagger(n,t_0)\eta$ source apply (`:391`) feeds `srcL[l] += W_ell[l][n]*rho`
(`:395`); the SAME sink `kphi`$=K(n,t)\phi'$ (`:411`) feeds `PhiL[l] += W_ell[l][n]*kphi` (`:414`).  Tower solves
`op_Dmsq.solve(psi_yl[IYL(l,b)], op_Dm(srcL[l]))` (`:400-402`); pair `Gyl[IYL(l,b)][dt] += psi_yl.dag(PhiL[l])`
(`:417`).  $(--)$ as in tp.  **Sign:** $\ell=0\approx0$, $\ell=1\sim-e^{-2t}$, $\ell=2\sim-e^{-3t}$ (Eqs. 4.34-4.35).

### 4.3 Spatial sp  (++ `:469-497`, parity (--) `:498-524`)

Insertion-diagonal over spatial links $a$ (`base.links[a]`, $il=$`map2il`):
$$
G^{s,++}(t_0\!\to\! t)=\frac1{4\pi}\sum_a w^{sp}_{il}\,\psi^{sp}_{a}(t_0)^\dagger\,K(\ell_a,t)\,\phi',\qquad
\psi^{sp}_a(t_0)=D_m^{-\dagger}K^\dagger(\ell_a,t_0)\eta .
$$
Code mirrors tp with `set_spatial(U,.,lk,.)` and `w_sp[il]`.  Own sink pass ($K(\ell_a,t)\phi'$, not shared with
the temporal pass).  $(--)$: tilde mirror (reuses `phimm`).  **Sign:** $\mathrm{Re}\,G_s>0$ (Eq. 4.28),
opposite to tp -- the cross-check.

---

## 5. Axial $C_{A+-}$  (code `:526-587`)

Vector + axial differ: own GW factor and $K^\dagger$ sink.  Valence legs (Sec. 1.1): flavor $m_F$ ->
massless $D_\text{ov}$ both legs; parity $m_P$ -> sink $\tilde D_{m_P}$; else $D_m$.  Forward leg and GW factor:
$$
\phi'_A=X^{-1}\eta\ (X=D_\text{ov}\text{ if flavor, else }D_m),\qquad \chi=(1-D_\text{ov})\,\phi'_A .
$$
Estimator (only $C_{A+-}$; $C_{A-+}=$ reflection $\Delta t\to N_t-\Delta t$, Eq. 3.57):
$$
C_{A+-}(t_0\!\to\! t)=\frac1{4\pi}\sum_a w_a\,\psi^A_a(t_0)^\dagger\,K^\dagger(a,t)\,\chi,\qquad
\psi^A_a(t_0)=X_\text{sink}^{-1}\,(1-D_\text{ov})\,K^\dagger(a,t_0)\,\eta,
$$
$X_\text{sink}=D_\text{ov}$ (flavor) / $\tilde D_{m_P}$ (parity) / $D_m$ (else).  Code: forward leg `:533-535`;
$\chi$ `:536`; tp source `:539-549` (`op_K(rho,eta)` -> `op_oneMinusD(rho,rho)` -> case solve into
`psi_tp[ITP(n,b)]`); tp sink `:551-558` (`set_temporal(U,t,n,dag=true)` -> `op_K(kphi,chi)` -> `Atp[b][dt] +=
w_tp[n]*psi_tp.dag(kphi)`); sp source/sink `:561-585` (same with links).  `psi_tp`/`psi_sp` are **reused**
(the vector results are already written).

---

## 6. Output  (h5; `corr.<k>.h5` in `data_<ESNID>/corr_nt0<N>_nhits<H>/`)

`ESNID = (free|<ens-basename>) + _vmRe<re>vmIm<im>`.  One file per config; resume skips iff the
`"complete"` sentinel (written LAST, `:474`-area) is present (read-only check `:367-376`).

Top-level: `t0s`, `n_t0`, `nhits`, `ls`, `complete`.  Per hit/origin:
```
h{h}/t0_{b}/tp/Vpp/{real,imag}        h{h}/t0_{b}/tp/Vmm/{real,imag}
h{h}/t0_{b}/sp/Vpp/{real,imag}        h{h}/t0_{b}/sp/Vmm/{real,imag}
h{h}/t0_{b}/ylm/Vpp/l{l}/{real,imag}  h{h}/t0_{b}/ylm/Vmm/l{l}/{real,imag}     l=0,1,2
h{h}/t0_{b}/axial/tp/Apm/{real,imag}  h{h}/t0_{b}/axial/sp/Apm/{real,imag}
```
Origins $t_0=b\cdot(N_t/n_{t_0})$, $b=0..n_{t_0}{-}1$; full $\Delta t=0..N_t{-}1$, index $(t-t_0)\bmod N_t$.
Downstream: average over $(h,t_0)$ per config, then jackknife configs.

---

## 7. Loop structure & solve count (per hit)

1. $\phi'=D_m^{-1}\eta$ -- 1 solve, shared.
2. Vector temporal: source solves $n_{t_0}(n_\text{sites}\text{[tp]}+n_\ell\text{[ylm]})$; **one** shared sink
   pass $K(n,t)\phi'$ ($n_\text{sites}N_t$ applies, fed to tp+ylm).
3. Vector spatial: source solves $n_{t_0}n_\text{links}$; sink pass $K(\ell,t)\phi'$.
4. Axial: forward leg + source solves $n_{t_0}(n_\text{sites}+n_\text{links})$; $K^\dagger$ sink passes.
5. Parity doubles the relevant source legs and adds a tilde sink pass.

Vector source solves/hit $=n_{t_0}(n_\text{sites}+n_\ell+n_\text{links})$ (identical to the three standalone
`jj_conn_*` programs); the unification computes the temporal sink applies **once** (shared tp+ylm, deduped over
$t_0$) instead of $2\,n_{t_0}\times$.
