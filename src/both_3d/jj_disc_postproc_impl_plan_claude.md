# Disconnected two-point post-processor -- implementation plan

## Goal

The measured current-current correlator is (PDF v2-10 Eqs. 3.39-3.40 -- note the **opposite**
relative sign between the two diagrams):
$$
\frac{2}{N_f}\langle J_{V,\pm}^{P}(0)\,J_{V,\pm}^{P}(t)\rangle_F
   \;=\; -\,C^P_{\pm\pm,c}(\Delta t) \;+\; C^P_{\pm\pm,d}(\Delta t),
\qquad P \in \{\,\text{tp},\ \text{sp},\ \text{ylm}\,\ell\,\}.
$$
The connected diagram (stored `Vpp`/`Vmm`) is $C_{++,c}=\mathrm{tr}(D_m^{-1}K(0)D_m^{-1}K(t))$
(Eq. 3.64); the disconnected diagram is the **product of single-current traces**
$C_{++,d}=\mathrm{tr}(D_m^{-1}K(0))\,\mathrm{tr}(D_m^{-1}K(t))=J(0)J(t)$ (Eq. 3.65), plain positive
product.  The full physical correlator therefore combines them as
$$
C^P_\text{full} \;=\; -\,C^P_\text{conn} \;+\; C^P_\text{disc}.
$$
**Convention (RESOLVED).** The stored `Vpp` $=C_{++,c}$ (Eq. 3.64), which asymptotes to a
**negative** curve in the data.  The physical correlator applies the explicit Eq. 3.39 minus:
$-\texttt{Vpp}$ ($>0$) is the connected contribution.  So the post-processor combines as
$-\,\texttt{conn}+\texttt{disc}$ (apply the minus to the connected leaf).  Disc enters with $+$.
Normalization stays at the bare $C$ level (no $2/N_f$; applied downstream).
`jj_corr_claude.cu` already writes, per config file `corr.<k>.h5`:

- **connected**, per source hit $h$ and origin $t_0$, indexed by $\Delta t$:
  `h{h}/t0_{b}/tp/Vpp`, `.../tp/Vmm`, `.../sp/{Vpp,Vmm}`, `.../ylm/{Vpp,Vmm}/l{l}` (length $N_t$).
- **disconnected single-current traces** (raw, origin-independent, length $N_t$, abs. time $t$):
  `h{h}/disc/tp/J`, `disc/sp/J`, `disc/ylm/l{l}/J`, and for parity `.../Jtil`.

This post-processor forms **only** the disconnected two-point $C^P_\text{disc}(\Delta t)$ from the
stored single-current traces $J^P$ (it does **not** read the connected keys), rearranging the raw
$J^P(t)$ into the projected, $\Delta t$-indexed two-point form that mirrors the connected leaves.
The $-C_\text{conn}+C_\text{disc}$ **sum is done downstream in the notebook** -- the connected
`corr.<k>.h5` is never read or mutated.  Done in a separate `.cc` (not the notebook): the disc
two-point is a per-config bilinear in $J$ with hit/translation bookkeeping -- cleaner compiled.

## Disconnected two-point estimator

The disconnected contribution to $\langle J^P(t_2) J^P(t_1)\rangle$ is the product of two
single-current loops on the **same gauge config**:
$$
C^P_\text{disc}(\Delta t) \;=\; \frac{1}{N_t}\sum_{t_0=0}^{N_t-1}
   \big\langle\, J^P(t_0+\Delta t)\; J^P(t_0) \,\big\rangle ,
$$
translation-averaged over all $t_0$ (the $J^P$ are full-$N_t$ and origin-independent, so all
$N_t$ origins are free statistics).

$J^P$ is a stochastic estimator over $n_h$ hits.  **Chosen estimator (decision 1): hit-averaged
product.**  Form $\bar J^P(t)=\frac1{n_h}\sum_h J^P_h(t)$, then
$$
C^P_\text{disc}(\Delta t) \;=\; \frac{1}{N_t}\sum_{t_0=0}^{N_t-1}\bar J^P(t_0+\Delta t)\,\bar J^P(t_0).
$$
This carries an $O(1/n_h)$ self-contraction bias (the loop variance), acceptable at large $n_h$.
The **unbiased** alternative (one-line switch) drops the $h=h'$ diagonal:
$\widehat{JJ}=\big[(\sum_h J_h(t_2))(\sum_h J_h(t_1))-\sum_h J_h(t_2)J_h(t_1)\big]/[n_h(n_h-1)]$
(needs $n_h\ge2$).

## Channel pairing (Vpp / Vmm)

Both current insertions are the **same** operator, so the disc loops pair like-with-like:

| channel | valence | connected legs | disc pairing |
|---|---|---|---|
| Vpp | massless / $m_F$ / $m_P$ | $D_m$ both legs | $J\!\cdot\!J$ |
| Vmm | massless / $m_F$ | (conj of Vpp) | $J^\ast\!\cdot\!J^\ast$ |
| Vmm | parity $m_P$ | tilde legs | $\tilde J\!\cdot\!\tilde J$ (`Jtil`) |

So: `Vpp` disc from `disc/*/J`; `Vmm` disc from `conj(J)` (non-parity) or `disc/*/Jtil` (parity).

## Output

Written **into the existing** `data_<ESNID>/corr_nt0<N>_nhits<H>/` dir as
`disc_proced.<k>.h5`, alongside the `corr.<k>.h5` it reads (so `corr.<k>.h5` is untouched).
**One length-$N_t$ array per config**, $\Delta t$-indexed, under the **same leaf names** as the
connected file so the notebook reuses its loader (just a different file prefix):
```
tp/Vpp   tp/Vmm    sp/Vpp   sp/Vmm    ylm/Vpp/l{l}   ylm/Vmm/l{l}        (= C^P_disc, length Nt, dt-indexed)
```
- `Vpp` = $\bar J(0)\bar J(t)$ from `disc/*/J`.
- `Vmm` = $\bar J^\ast(0)\bar J^\ast(t)$ (non-parity, $=\mathrm{conj}(Vpp)$) or
  $\bar{\tilde J}(0)\bar{\tilde J}(t)$ from `disc/*/Jtil` (parity).
- No axial (disc has no axial loop -- axial is connected-only).
- `complete` sentinel written last (same resume gate; a separate sentinel from `corr.<k>.h5`'s).

Downstream the notebook forms the physical correlator per channel/leaf as
$-\,\texttt{corr.Vpp} + \texttt{disc\_proced.Vpp}$.

## Files

- NEW `jj_disc_postproc_claude.cc` -- the post-processor (HDF5 read/write via HighFive, as in the `.cu`).
- EDIT `run_jj_analysis_claude.sh` -- add a post-proc step (loop configs in each `corr_nt0...` dir after `--all`).
- EDIT `jj_conn_tpproj_various_mass_claude.ipynb` -- loader for `disc_proced.<k>.h5` + the $-$conn$+$disc sum cell.
- This plan: `jj_disc_postproc_impl_plan_claude.md`.

## Ordered chunks

1. **Reader.** Iterate `corr.*.h5` in a `corr_nt0...` dir; per config read the disc traces `h{h}/disc/{tp,sp,ylm/l{l}}/{J,Jtil}` (skip `h{h}/t0_*` connected keys entirely); detect parity from presence of `Jtil`.  Files: `jj_disc_postproc_claude.cc`.
2. **Disc two-point.** Hit-average $\bar J^P$ over $h$; translation-average $\bar J^P(t_0+\Delta t)\bar J^P(t_0)$ over all $t_0$ -> $C^P_\text{disc}(\Delta t)$, both channels (Vpp from `J`, Vmm from `conj(J)` or `Jtil`).  Files: same.
3. **Write.** `disc_proced.<k>.h5` in the same dir, leaves `tp/{Vpp,Vmm}`, `sp/{Vpp,Vmm}`, `ylm/{Vpp,Vmm}/l{l}`; `complete` sentinel; resume-skip on it.  Files: same.
4. **Driver + notebook.** Run-script post-proc step; notebook `disc_proced` loader + the $-$conn$+$disc plot cell.  Files: `run_jj_analysis_claude.sh`, notebook.

## Resolved decisions

1. **Hit-bias:** hit-averaged product $\bar J(t_2)\bar J(t_1)$, $\bar J=\frac1{n_h}\sum_h J_h$ (one-line switch to cross-hit if wanted later).
2. **Translation average:** over all $N_t$ origins.
3. **Output granularity:** one $(h,t_0)$-averaged array per config; notebook jackknifes configs.
4. **Sign / norm (Q6):** $C_\text{full}=-\,C_\text{conn}+C_\text{disc}$ at the bare $C$ level (no $2/N_f$); stored `Vpp`$=C_{++,c}$.  Sum done in the notebook.
5. **Output location (A):** `disc_proced.<k>.h5` **inside** the existing `corr_nt0<N>_nhits<H>/` dir; holds the **processed disc two-point only** (no conn read, no in-file sum).  `corr.<k>.h5` untouched.

All decisions resolved -- ready to code on go-ahead.
