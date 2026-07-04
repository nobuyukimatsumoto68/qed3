# Audit -- scalar $Y_{\ell m}$ correlator drivers (connected + disc)

Audit of the implementation of `scalar_ylm_corr_impl_plan_claude.md` in:

- `jj_local_ylm_scalar_conn_stoch_claude.cu` (connected $V_{++}$ tower)
- `jj_local_ylm_scalar_disc_stoch_claude.cu` (disc one-point loop $J^0_{\ell m}$ + condensate)

Both are copies of the vector/axial drivers (`jj_local_ylm_conn_stoch_claude.cu`,
`jj_local_ylm_disc_stoch_claude.cu`) with an $a=0$ identity-vertex channel and `--IsScalarOnly` added.
Focus per user request: correctness, and readability of the conventions inherited from the vector code
(e.g. `Vpp`/`Vmm`) that make the scalar path hard to debug.

## Verdict

Physics is correct in both files. The estimators match the plan derivation and reduce to the intended
traces. No hard correctness bug found. The issues are **convention/normalization asymmetries between the two
drivers and versus the plan text**, which are latent traps for the (still-to-be-written) analysis notebook,
plus a few readability warts.

## Physics -- verified correct

### Connected (`jj_local_ylm_scalar_conn_stoch_claude.cu`)

- `phi = D_m^{-1} eta`: `op_DmH` then `op_Dmsq.solve` $=(D^\dagger D)^{-1}D^\dagger\eta$ (`:357-358`). OK
- `psi = D_m^{-dag} src`: `op_Dmsq.solve` then `op_Dm` $=(D^\dagger)^{-1}W_{\ell m}(t_0)\eta$ (`:400-401`). OK
- `gvs[dt] += psi.slice(t).dot(Phi.slice(t))` $=\psi^\dagger W_t\phi=\eta^\dagger W_0 D^{-1}W_t D^{-1}\eta$, which
  under $E[\eta\eta^\dagger]=\mathbb 1$ estimates $\operatorname{tr}[W_0 D^{-1}W_t D^{-1}]=V_{++,\ell m}(t)$
  (`:406-408`). OK -- matches plan.
- Identity vertex (no `mult_sigma`), both $W$'s carry `mult_Ylm_real` (folds $A_n Y_{\ell m}$), same $(\ell,m)$
  at source and sink; spin+space traced together -- correct for the scalar bilinear (`:398,:404`). OK
- Vector+axial block is byte-identical to the parent (verified by `diff`); no regression.

### Disc (`jj_local_ylm_scalar_disc_stoch_claude.cu`)

- Scalar loop `JS[l][m] += accumulate_loop_raw(W_lm phi)` with `Gamma_phi = mult_Ylm_real(phi)` (identity vertex,
  no `mult_sigma`) estimates $J^0_{\ell m}(t)=\operatorname{tr}[W_{\ell m}(t)D_m^{-1}]$, stored COMPLEX
  (`:369-376,:413-419`). Analysis takes $\mathrm{Re}\to\sigma_{PS}$, $\mathrm{Im}\to\sigma_{FS}$. OK -- matches plan.
- Spin loop runs `spin=0,1` regardless of channel, so the $a=0$ loop is traced over both spin components. OK
- The `assert(!parity)` correctly forbids the out-of-scope $m_P$ leg in both files (conn `:227-228`, disc `:213-214`).

## Findings (ranked)

### F1 [HIGH -- analysis trap] conn is $1/(4\pi)$-folded, disc is raw, and the plan text says the opposite

- Connected `write_corr` multiplies every stored trace by `inv4pi` at MEASUREMENT time
  (conn `:295-300`, "matches jj_local_deter").
- Disc `write_vec` stores $J^0$ RAW, no `inv4pi` (disc `:282-287`, "matching disc_claude.cu").
- The plan states the "$1/(4\pi)$ fold ... applied in analysis (measurement stores raw per-config traces)"
  (plan `:113-114`) -- which is true for the disc but FALSE for the connected.

Consequence: the analysis reconstructing $\frac{N_f}{2}(-C_c+C_d)$ must NOT re-fold $C_c$ (already folded) but
MUST fold $C_d$ (raw) so the two terms share one $Y_{\ell m}$ normalization. This asymmetry is inherited from
the validated vector conn/disc pair, so it is presumably "correct if you mirror the vector analysis" -- but a
fresh scalar notebook following the plan literally will either double-fold $C_c$ ($\times 1/4\pi$ error) or
leave $C_d$ unfolded (relative $4\pi$ between connected and disc). Fix: correct the plan text and add a one-line
"already inv4pi-folded" / "RAW" banner at each writer so the two conventions are explicit.

### F2 [MED -- layout] the two drivers disagree on where the scalar tower lives

- Connected scalar path: `h0/scalar/l{l}/m{m}/{Vpp,Vmm}` -- NOT under `ylm/`, channel named `scalar`
  (conn `:453`).
- Disc scalar path: `h0/disc/ylm/s0/l{l}/m{m}/J` -- under `ylm/`, channel named `s0`
  (disc `:416`), consistent with its own vector `h0/disc/ylm/s{a}/`.

So the same $a=0$ channel is `scalar` in one file and `s0` in the other, and lives under different parent
groups. Analysis must special-case both. Recommend one scheme, e.g. connected `h0/ylm/s0/l/m/` (mirrors its
vector `h0/ylm/s{a}/`) or disc `h0/disc/scalar/` -- but pick one $a=0$ name and use it in both.

### F3 [MED -- readability, the `Vpp`/`Vmm` wart you flagged] connected stores a redundant conjugate

`write_corr(kp+"Vpp", gS, false)` and `write_corr(kp+"Vmm", gS, true)` write the SAME accumulator `gS`, the
second as its complex conjugate (conn `:454-455`). So `Vmm` $\equiv \operatorname{conj}(V_{++})$ -- zero new
information, double the scalar HDF5 output. In the vector/parity context `Vpp`/`Vmm` were genuinely distinct
$V_{++}/V_{--}$ channels; here the name misleads a debugger into thinking two correlators exist.

Note the DISC driver already does the clean thing: it stores ONE complex $J^0$ and lets analysis take
Re/Im. The connected should adopt the same pattern -- store a single complex `V` (or `Vpp`) and form
$2\,\mathrm{Re}$ in analysis -- rather than persisting a conjugate copy under a parity-channel name.

### F4 [MED -- sign documentation] `C_c`/`-C_c` conflated; loop-minus asymmetry between $J^0$ and the condensate

- Conn comment `:450` reads `C_c = -(Vpp+Vmm)`, but the plan defines $C_c=2\,\mathrm{Re}\,V_{++}=\texttt{Vpp}+\texttt{Vmm}$
  (positive); the minus belongs to the full $\frac{N_f}{2}(-C_c+C_d)$, not to $C_c$. Stored `Vpp` is $+V_{++}$.
- In the disc file, the condensate carries the closed-fermion-loop minus: `etadag_xi = -acc_etadag_xi`
  (disc `:410`), whereas the $J^0$ tower is stored RAW with NO minus (disc `:374`, `accumulate_loop_raw`).
  This is actually CORRECT for the two-point $C_d=\texttt{two\_point}(J^0)$ because the two loop minus signs
  square away -- but a reader comparing the $\ell=0$ $J^0$ tower against `condensate/etadag_xi` will find them
  off by a sign and chase a phantom bug. Both facts deserve an explicit comment.

### F5 [LOW -- CLI] `--IsScalarOnly` is CamelCase

Every other flag is hyphen-lowercase (`--spin-dilution`, `--mass-re`, `--disc-tblock`); `--scalar-only`
would match (conn `:177`, disc `:165`).

### F6 [LOW -- accuracy caveat] Zolotarev pole count fixed at 11

`Fermion Dm(DW, valence_mass, 11)` with no fixed window (conn `:268`, disc `:253`), same as the parent
drivers. Fine at $L=1$; the project's $L\ge 2$ work needed npole 21 + window 0.001 near zero modes. If either
driver is run at `N_REFINE_CLI >= 2`, valence $D^{-1}$ accuracy is a caveat worth a comment.

### F7 [LOW -- stale comment] connected work-vector comment

`FermionVector ... // psi = vector source leg; psiA = axial source leg` (conn `:302`) does not mention that
`psi` is reused as the scalar source leg in the `:399` block.

## Not a bug -- verified

- `--IsScalarOnly` append path in both files: copies the complete file, opens ReadWrite, adds only the scalar
  keys (no key collision with `h0/ylm/...` or `h0/disc/...`), skips re-writing `complete` to avoid a dup-key
  crash, and gates reruns on scalar-key presence (conn `:325-337,:420-460`; disc `:307-323,:384-420`).
  The disc `has_s0` check probes a concrete dataset (`h0/disc/ylm/s0/l0/m0/J/real`, disc `:315`) -- more
  robust than the conn `getGroup("h0").exist("scalar")` (conn `:330`).
- `.tmp` name is per-$(k,h)$, so it is safe under the project's disjoint-config MPS packing (no two processes
  write the same $(k,h)$).
- Scalar block runs in default (full) mode too, adding the scalar tower alongside vector+axial in one file --
  intended per plan.

## Recommended edits (in priority order)

1. F1: correct plan `:113-114`; add "already inv4pi-folded" / "RAW -- fold in analysis" banners at conn
   `write_corr` and disc `write_vec`.
2. F2: unify the $a=0$ HDF5 path/name across the two drivers.
3. F3: drop `Vmm` in the connected scalar output; store one complex `V`, do $2\,\mathrm{Re}$ in analysis.
4. F4: fix the `C_c = -(Vpp+Vmm)` comment; document the $J^0$-vs-condensate loop-sign convention.
5. F5-F7: rename flag, add npole caveat, refresh the work-vector comment.
