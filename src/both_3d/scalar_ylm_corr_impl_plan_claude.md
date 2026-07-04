# Impl plan -- scalar-density $Y_{\ell m}$ correlators $\sigma_{FS},\sigma_{PS}$ (connected + disc reconstruction)

Add the spherical-harmonic tower of the SCALAR-density two-point functions $\langle\sigma_{PS}\sigma_{PS}\rangle$
and $\langle\sigma_{FS}\sigma_{FS}\rangle$, using the SAME $Y_{\ell m}$ machinery as the local vector/axial
current tower. Copy the existing production ylm driver and gate the scalar-only mode with `--IsScalarOnly`.

Reference note: `qed3int_v2-14.pdf` (source `qed3int_v2(7)/main.tex`). Copy targets:
`jj_local_ylm_conn_stoch_claude.cu` (connected), `jj_local_ylm_disc_stoch_claude.cu` (disc one-point loops).
Estimator/Ylm machinery: `jj_local_ylm_impl_plan_claude.md`. Densities audited in
`jj_corr_dilute_audit_claude.md` Sec. 2 and `condensate_contact_massive_claude.md`.

## IMPORTANT -- the scalar CORRELATORS are NOT in the note pdf

The PDF defines only the scalar **densities** (mass-term bilinears):
$$
\sigma_{PS} = \eta^\dagger\xi + \xi^\dagger\eta = \bar\Psi\Psi \qquad(\text{Eq. 1.23; couples to } m_F\in\mathbb R,\ \text{Eq. 2.1})
$$
$$
\sigma_{FS} = \eta^\dagger\xi - \xi^\dagger(1-D_{\rm ov}^\dagger)\eta \qquad(\text{Eq. 1.55; couples to } m_P\in i\mathbb R,\ \text{Eq. 2.2})
$$
$\sigma_{PS}$ = parity-symmetric / flavor-breaking; $\sigma_{FS}$ = flavor-symmetric / parity-breaking; the
Ginsparg-Wilson unitary is $V=D_{\rm ov}-1$, so $\xi^\dagger V^\dagger\eta = -\xi^\dagger(1-D_{\rm ov}^\dagger)\eta$.
The **two-point** correlators, their connected/disc Wick structure, and any $Y_{\ell m}$/sp/tp projection or
scalar weight are **absent from v2-14** (verified: only the current $C_V/C_A$ tower, Eqs. 3.39-3.66, 4.27-4.36,
is written out; "condensate" is a one-line idea in the Discussion). The Wick derivation below is therefore
**mine** and must be confirmed before coding.

## Propagator legs -- Eqs. (3.57)/(3.58), massless / $m_F$ only ($m_P$ dropped, per user)

Nonzero legs only ($\xi,\eta$ couple only off-diagonally, so $\langle\xi\xi^\dagger\rangle=\langle\eta\eta^\dagger\rangle=0$):
$$
P \equiv \langle\xi\,\eta^\dagger\rangle = D_m^{-1}\ (\text{Eq. 3.57, forward}),\qquad
\bar P \equiv \langle\eta\,\xi^\dagger\rangle = (D_m^{-1})^\dagger = D_m^{-\dagger}\ (\text{Eq. 3.58, backward}),\qquad D_m=D_{\rm ov}+m .
$$

## Wick derivation -- CONNECTED **and** DISCONNECTED (transcribed from `scalar_handwritten.pdf`, user, 2026-07-03)

**Correction to my earlier claim.** I had said $\langle T_1 T_2\rangle=0$. That is only true for its
*connected* contraction; the $T_1 T_2$ product IS nonzero as a **disconnected** loop product (each of $T_1,T_2$
self-contracts into a one-point trace). So the full scalar two-point has BOTH a connected and a disconnected
term, with an overall flavor factor $N_f/2$. My connected pieces below were right; the disconnected piece was
missing.

**FURNISHED definition (user, 2026-07-03 -- REVERTED the earlier $(1-D_{\rm ov})\to1$): both insertions carry
$(1-D_{\rm ov})$, EXACTLY as the handwritten note.** $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{\rm ov}^\dagger)\eta$
(the full lattice Eq. 1.55). The wall source at $t_0$ still works (localize with $W_{\ell m}(t_0)$ on $\eta$
FIRST, then apply $(1-D_{\rm ov})$, then solve); the $(1-D_{\rm ov})$ ride the operator chain, not the source
localization. Cost: the FS channel adds one shared sink-base inversion $\chi$ + a per-$(\ell,m)$ source solve
$\psi^{FS}$ (roughly 2x the connected solves for FS); $\sigma_{PS}$ is UNCHANGED (never had a GW factor).
$Y_{\ell m}$-leg convention: the harmonic weight attaches to the $\eta$-leg point (where $D_m^{-\dagger}$ meets
$(1-D_{\rm ov}^\dagger)$), held CONSISTENT between conn and disc; the other-leg choice differs by $O(a)$.

Scalar vertex is the **identity** ($a=0$, no Pauli $\sigma_a$); harmonic-weighted insertion
$W_{\ell m}(t)=\sum_n A_n Y_{\ell m}(\hat n)\,\mathbb 1_{\rm spin}(n,t)$, $A_n=\texttt{dual\_areas}[n]$.
Write $\sigma_{PS}=T_1+T_2^{PS}$, $\sigma_{FS}=T_1+T_2^{FS}$: $T_1=\eta^\dagger\xi$, $T_2^{PS}=\xi^\dagger\eta$,
$T_2^{FS}=-\xi^\dagger(1-D_{\rm ov}^\dagger)\eta$ (FURNISHED; same two terms as the one-point condensate
`etadag_xi`, `xidag_1mDdag_eta`).

### $\sigma_{PS}$ (page 1 left)
$$
\big\langle(\eta^\dagger\xi+\xi^\dagger\eta)(\eta^\dagger\xi+\xi^\dagger\eta)'\big\rangle
=\underbrace{\langle(\eta^\dagger\xi)(\eta^\dagger\xi)'\rangle_G+\langle(\xi^\dagger\eta)(\xi^\dagger\eta)'\rangle_G}_{\text{connected}}
+\underbrace{\big\langle(\eta^\dagger\xi+\xi^\dagger\eta)(\eta^\dagger\xi+\xi^\dagger\eta)'\big\rangle_G^{\rm loops}}_{\text{disconnected}}
$$
$$
=\frac{N_f}{2}\Big\langle
-\operatorname{tr}[D_m^{-1}(x,x')D_m^{-1}(x',x)]-\operatorname{tr}[D_m^{-\dagger}(x,x')D_m^{-\dagger}(x',x)]
+\big(\operatorname{tr}D_m^{-1}(x,x)+\operatorname{tr}D_m^{-\dagger}(x,x)\big)\big(\operatorname{tr}D_m^{-1}(x',x')+\operatorname{tr}D_m^{-\dagger}(x',x')\big)\Big\rangle_U
$$
$$
=\frac{N_f}{2}\Big\langle
-2\,\mathrm{Re}\operatorname{tr}[D_m^{-1}(x,x')D_m^{-1}(x',x)]
+4\,[\mathrm{Re}\operatorname{tr}D_m^{-1}(x,x)]\,[\mathrm{Re}\operatorname{tr}D_m^{-1}(x',x')]\Big\rangle_U
\equiv\frac{N_f}{2}\big\{-\langle C_{PS,c}(U)\rangle_U+\langle C_{PS,d}(U)\rangle_U\big\}.
$$
$$
C_{PS,c}=2\,\mathrm{Re}\operatorname{tr}[D_m^{-1}(x,x')D_m^{-1}(x',x)],\qquad
C_{PS,d}=4\,[\mathrm{Re}\operatorname{tr}D_m^{-1}(x,x)]\,[\mathrm{Re}\operatorname{tr}D_m^{-1}(x',x')].
$$

### $\sigma_{FS}$ (page 1 right, FURNISHED: $(1-D_{\rm ov})$ on both insertions, per the note)
$$
\big\langle(\eta^\dagger\xi-\xi^\dagger(1-D_{\rm ov})\eta)(\cdots)'\big\rangle
=\underbrace{\langle(\eta^\dagger\xi)(\eta^\dagger\xi)'\rangle_G+\langle(\xi^\dagger(1-D_{\rm ov})\eta)(\xi^\dagger(1-D_{\rm ov})\eta)'\rangle_G}_{\text{connected}}
+\underbrace{\langle(\cdots)(\cdots)'\rangle_G^{\rm loops}}_{\text{disconnected}}
$$
$$
=\frac{N_f}{2}\Big\langle
-\operatorname{tr}[D_m^{-1}(x,x')D_m^{-1}(x',x)]
-\operatorname{tr}[(1-D_{\rm ov}^\dagger)D_m^{-\dagger}(x,x')\,(1-D_{\rm ov}^\dagger)D_m^{-\dagger}(x',x)]
$$
$$
\qquad
+\big[\operatorname{tr}D_m^{-1}(x,x)-\operatorname{tr}((1-D_{\rm ov})D_m^{-1}(x,x))\big]\big[\operatorname{tr}D_m^{-1}(x',x')-\operatorname{tr}((1-D_{\rm ov})D_m^{-1}(x',x'))\big]\Big\rangle_U
\equiv\frac{N_f}{2}\big\{-\langle C_{FS,c}\rangle+\langle C_{FS,d}\rangle\big\}.
$$
$T_1 T_1$ gives the shared $V_{++}=\operatorname{tr}[D_m^{-1}D_m^{-1}]$; $T_2^{FS}T_2^{FS}$ gives the GW-dressed
$V_{--}^{FS}=\operatorname{tr}[\Gamma\Gamma]$, $\Gamma=(1-D_{\rm ov}^\dagger)D_m^{-\dagger}$. Massless / $m_F$; the
$(1-D_{\rm ov})$ is the massless overlap. ($D_m^{-1}$ and $(1-D_{\rm ov})$ are both functions of $D_{\rm ov}$ so
they commute.)

### Stochastic estimators (one-end trick, $\phi=D_m^{-1}\eta$, $Z_2$ wall at $t_0$; per $\ell,m$; $W$ identity vertex)
**Connected $\sigma_{PS}$** ($a=0$ `Vpp`; also the $T_1T_1$ term shared by $\sigma_{FS}$):
$$
\psi_{\ell m}=D_m^{-\dagger}W_{\ell m}(t_0)\eta,\quad
g_{\ell m}(t)=\psi_{\ell m}^\dagger[W_{\ell m}(t)\phi]\to V_{++,\ell m}=\operatorname{tr}[W_0 D_m^{-1}W_t D_m^{-1}];\quad
C^{PS}_c=V_{++}+V_{++}^*.
$$
Store $V_{++}$ as `h0/scalar/l{l}/m{m}/{Vpp,Vmm=conj}`.
**Connected $\sigma_{FS}$ extra term** ($V_{--}^{FS}$, GW on both insertions):
$$
\Gamma=(1-D_{\rm ov}^\dagger)D_m^{-\dagger},\quad
\chi=\Gamma\eta=(1-D_{\rm ov}^\dagger)D_m^{-\dagger}\eta\ (\text{1 solve/hit, shared}),\quad
\psi^{FS}_{\ell m}=\Gamma^\dagger W_{\ell m}(t_0)\eta=D_m^{-1}(1-D_{\rm ov})W_{\ell m}(t_0)\eta,
$$
$$
g^{FS}_{\ell m}(t)=\big(\psi^{FS}_{\ell m}\big)^\dagger[W_{\ell m}(t)\chi]\to V_{--,\ell m}^{FS}=\operatorname{tr}[W_0\Gamma W_t\Gamma];\quad
C^{FS}_c=V_{++}+V_{--}^{FS}.
$$
Store $V_{--}^{FS}$ as `h0/scalar_fs/l{l}/m{m}/Vmm`.
**Disconnected** -- two COMPLEX one-point loops on the disc $\phi$ (NO extra solve; $(1-D_{\rm ov})\phi$ is a
cheap mat-vec): $J_1=\operatorname{tr}[W_{\ell m}(t)D_m^{-1}]$ and $J_{1mD}=\operatorname{tr}[W_{\ell m}(t)(1-D_{\rm ov})D_m^{-1}]$.
$$
J^{PS}_{\ell m}=J_1+J_1^{*}=2\,\mathrm{Re}\,J_1,\qquad
J^{FS}_{\ell m}=J_1-J_{1mD}^{*},\qquad
C_{d}\leftrightarrow\tfrac{1}{2\ell+1}\textstyle\sum_m\texttt{two\_point}(J^{\cdot}_{\ell m}).
$$
Store $J_1$ as `h0/disc/ylm/s0/.../J` and $J_{1mD}$ as `h0/disc/ylm/s0_1mD/.../J`.
Analysis: $g_\ell=\frac{1}{2\ell+1}\sum_m g_{\ell m}$; $1/(4\pi)$ fold + $N_f/2$ and $-C_c+C_d$ applied in analysis.

Cost per hit (per spin class): $\sigma_{PS}$ conn = $1(\phi)+16(\psi)$. $\sigma_{FS}$ conn ADDS $1(\chi)+16(\psi^{FS})$.
Disc: $J_1$ and $J_{1mD}$ ride the disc $\phi$ (one $(1-D_{\rm ov})\phi$ mat-vec, no extra solve).

## Disc piece -- one-point trace towers (Chunk 3, FURNISHED)

Two COMPLEX one-point loops on the disc $\phi=D_m^{-1}\eta$: $J_1=\operatorname{tr}[W_{\ell m}D_m^{-1}]$ (identity,
`h0/disc/ylm/s0`) and $J_{1mD}=\operatorname{tr}[W_{\ell m}(1-D_{\rm ov})D_m^{-1}]$ (GW, `h0/disc/ylm/s0_1mD`;
apply $(1-D_{\rm ov})$ to $\phi$ then the $W_{\ell m}$ readout -- a cheap mat-vec, NO extra solve). Analysis:
$J^{PS}=2\mathrm{Re}\,J_1$, $J^{FS}=J_1-J_{1mD}^{*}$; $C_d=\texttt{two\_point}(J)$.

## Files to modify / create

- NEW `jj_local_ylm_scalar_conn_stoch_claude.cu` (copy of `jj_local_ylm_conn_stoch_claude.cu`): add the $a=0$
  scalar connected channel (`Vpp`/`Vmm=conj`; ONE computation, shared by PS and FS after $(1-D_{\rm ov})\to1$),
  add `--IsScalarOnly`, add the append-mode writer + scalar-key completion gate. Vector/axial path kept.
- NEW `jj_local_ylm_scalar_disc_stoch_claude.cu` (COPY of `jj_local_ylm_disc_stoch_claude.cu`; original left
  UNTOUCHED): add the $a=0$ COMPLEX scalar tower $J^0_{\ell m}(t)$ -> `h0/disc/ylm/s0/l{l}/m{m}/J` +
  `--IsScalarOnly` append.
- NEW run script(s) `run_ylm_scalar_*_claude.sh`.
- Analysis: reconstruct $C^{PS/FS}=\frac{N_f}{2}(-C_c+C_d)$ per-$\ell$ (new/updated notebook, NOT overwritten).

## Ordered chunks

1. **Scalar connected driver** -- copy conn driver; add the $a=0$ connected channel (`Vpp`/`Vmm=conj`) via the
   existing estimator with `mult_sigma` skipped (identity vertex), only `mult_Ylm_real`. ONE channel serves
   both PS and FS. Output `h0/scalar/l{l}/m{m}/{Vpp,Vmm}`. Files: `jj_local_ylm_scalar_conn_stoch_claude.cu`.
2. **`--IsScalarOnly` + append mode** -- flag off: compute vector+axial (as now) AND scalar in one file;
   flag on: compute ONLY scalar and APPEND `h0/scalar/...` into the existing conn per-hit `.h5` (atomic
   copy-to-`.tmp`+rename; skip only if the scalar keys already exist, NOT the `complete` gate). Files: same driver.
3. **Disc scalar tower** -- COPY `jj_local_ylm_disc_stoch_claude.cu` -> `jj_local_ylm_scalar_disc_stoch_claude.cu`
   (original untouched), add the $a=0$ COMPLEX one-point loop $J^0_{\ell m}(t)$ -> `h0/disc/ylm/s0/l{l}/m{m}/J`
   (Re->PS, Im->FS in analysis) + `--IsScalarOnly` append. Files: `jj_local_ylm_scalar_disc_stoch_claude.cu`.
4. **Run scripts + analysis** -- reconstruct $-C_c+C_d$ per-$\ell$. Files: `run_ylm_scalar_*_claude.sh`, notebook.

## Decisions (user, 2026-07-03) -- IDENTITY-ERA, SUPERSEDED by "REVERT to FURNISHED" + "Current status" below

(Kept as history; the $(1-D_{\rm ov})\to1$ bullet is NO LONGER the plan -- see the FURNISHED derivation at top
and the Current status block at the very end.)

- ~~**Both** $\sigma_{PS}$ and $\sigma_{FS}$, with $(1-D_{\rm ov})\to1$ (continuum form, Eq. 1.29) so the vertex is
  LOCAL and wall-sourceable. Consequence: connected identical for PS/FS (one computation); disc differs
  (Re vs Im of the same one-point loop).~~  -- REVERTED to furnished (both insertions carry $(1-D_{\rm ov})$).
- Full correlator $=\frac{N_f}{2}(-C_c+C_d)$; $C_c=2\mathrm{Re}\operatorname{tr}[D_m^{-1}D_m^{-1}]$,
  $C_d=$ loop product ($+$Re for PS, $-$Im$^2$ for FS). $\langle T_1T_2\rangle$ enters the DISC (loop product),
  not the connected -- my earlier "=0" was wrong for the disconnected part.
- Legs = Eqs. **(3.57)/(3.58)**; **massless + $m_F$ (real) only**, NO parity ($m_P$/tilde-$D$).
- `--IsScalarOnly`: **APPEND** `h0/scalar/...` into the existing conn per-hit `.h5` (atomic copy-to-`.tmp`+rename);
  default (flag off) computes vector+axial AND scalar in one file.
- Disc: **COPY** `jj_local_ylm_disc_stoch` -> `jj_local_ylm_scalar_disc_stoch_claude.cu` (original untouched,
  house rule) and add the $a=0$ complex $J^0_{\ell m}(t)$ tower there.
- Scalar vertex is the **identity** (no Pauli $\sigma_a$).

## REVERT to FURNISHED (user, 2026-07-03) -- supersedes the $(1-D_{\rm ov})\to1$ version

Bring back the full GW factor on BOTH insertions (Eq. 1.55). Sub-tasks on top of the shipped identity code:
- R2 conn: add massless $D_{\rm ov}$ apply; compute $\chi=(1-D_{\rm ov}^\dagger)D_m^{-\dagger}\eta$ (1 solve/hit) and
  per-$(\ell,m)$ $\psi^{FS}=D_m^{-1}(1-D_{\rm ov})W_{\ell m}(t_0)\eta$; new output `h0/scalar_fs/l/m/Vmm`
  ($V_{--}^{FS}$). `h0/scalar` (Vpp/Vmm=conj) stays = $\sigma_{PS}$. `--IsScalarOnly` append gains scalar_fs keys.
- R3 disc: add the second loop $J_{1mD}=\operatorname{tr}[W(1-D_{\rm ov})D_m^{-1}]$ -> `h0/disc/ylm/s0_1mD`
  ((1-D_ov)phi mat-vec, no extra solve).
- R4 notebook: $\sigma_{PS}$ unchanged; $\sigma_{FS}$ conn $=-(N_f/2)(V_{++}+V_{--}^{FS})$ (Vpp + scalar_fs/Vmm),
  disc from $J_1-J_{1mD}^{*}$.

## Progress (2026-07-03) -- IDENTITY-ERA chunks (SUPERSEDED; see "Current status" at the very end)

- Chunk 1 (conn scalar $V_{++}$/`Vmm=conj`): DONE, smoke PASS.
- Chunk 2 (conn `--IsScalarOnly` + append): DONE, smoke PASS (Phase B append verified).
- Chunk 3 (disc scalar $J^0$ tower, NEW copy file): DONE, smoke PASS; original disc driver reverted (git-clean).
- Chunk 4 (run script + analysis): DONE.
  - `run_ylm_scalar_claude.sh [together|scalaronly]` -- env knobs GSQ/NREF/KMIN/STRIDE/NCONF/DISC_TB;
    `together`=fresh vector+axial+scalar (new ensembles), `scalaronly`=`--IsScalarOnly` append onto existing
    jj files.  DISC_TB default 2 (must match existing disc dir for append).
  - `jj_scalar_ylm_analysis_claude.ipynb` (NEW, first creation) -- reconstructs $\tfrac{N_f}{2}(-C_c+C_d)$;
    $C_c=2\mathrm{Re}\,V_{++}$ from `h0/scalar`, $C_d$ from `h0/disc/ylm/s0` ($2\mathrm{Re}\,J^0\to$PS,
    $2i\,\mathrm{Im}\,J^0\to$FS); mirrors the vector notebook (jackknife, `two_point_pp`, Eq.12 DC-sub, plateau-sub).

All 4 chunks implemented and PASS (Chunks 1-4, user-confirmed).  Normalization audited OK (conn folds
inv4pi / disc raw; nb re-applies inv4pi ONCE).  Remaining: production via `run_ylm_scalar_claude.sh` on the
interacting ensembles, then the analysis notebook on real data.

## Current status (FURNISHED, 2026-07-03) -- THE authoritative status

The plan is FURNISHED: $(1-D_{\rm ov})$ on both $\sigma_{FS}$ insertions (Eq. 1.55). All identity-era text above is
history. Implemented + consistent across `.md` -> drivers -> smoke -> notebook -> memory:
- CONN `jj_local_ylm_scalar_conn_stoch_claude.cu`: `h0/scalar/{Vpp,Vmm=conj}` ($\sigma_{PS}$) +
  `h0/scalar_fs/Vmm` ($V_{--}^{FS}$; separate massless $D_{\rm ov}$, shared $\chi$ + per-mode $\psi^{FS}$).
- DISC `jj_local_ylm_scalar_disc_stoch_claude.cu`: `h0/disc/ylm/s0` ($J_1$) + `h0/disc/ylm/s0_1mD` ($J_{1mD}$).
- Notebook `jj_scalar_ylm_analysis_claude.ipynb` (9 cells): `Cc_PS=2\mathrm{Re}V_{++}`, `Cc_FS=V_{++}+V_{--}^{FS}`;
  disc PS$=2\mathrm{Re}J_1$, FS$=J_1-\mathrm{conj}(J_{1mD})$.
- Smoke scripts updated (masses 0.07/0.09; check `scalar_fs`/`s0_1mD`).
- INDEPENDENT AUDIT of furnished code vs this `.md` = **CLEAN** (from-scratch Wick re-derivation; confirmed the
  code builds ordered $\Gamma$ explicitly and does NOT rely on $[D_m^{-1},(1-D_{\rm ov})]=0$; separate massless
  $D_{\rm ov}$). 3 LOW notes (user said leave for now): (i) conn `--IsScalarOnly` gate tests `h0/scalar` not
  `scalar_fs` (only bites pre-revert legacy files); (ii) `op_Dmsq` DDH normality (pre-existing, exact at m=0);
  (iii) minor doc shorthand.
- Both original drivers verified git-clean.

NEXT: re-run the two smoke scripts (furnished keys), then production via `run_ylm_scalar_claude.sh`.
