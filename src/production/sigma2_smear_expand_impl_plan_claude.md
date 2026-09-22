# $N_v$-truncation (distillation-smearing) GEVP expansion of the PS $\sigma^2$ point basis -- impl plan

## Goal / physics
Expand the established PS$^2$ $\{\sigma^2_{00}, O_{2m}, O_{1m}\}$ Hankel+rebase GEVP (Nf2 gsq1.0 L1, 400 cfg _v1,
state0$\approx0.46$ one-meson-rich / state1$\approx0.62$ two-meson) with a **multi-smearing variational basis**
built by TRUNCATING the distillation perambulator to $N_v\in\{4,8,16\}$ modes (plus 24 = full = the existing
anchor). More interpolators with different spatial profiles -> better resolution of the $0^{++}$ tower and,
hopefully, a sharper two-meson $m_1$ (the sub-threshold / bound-state question). Later: the flavor-cross matrix
$\{\sigma_{PS}^2,\sigma_{FS}^2,\sigma_{FS}\sigma_{PS}\}$ where FS$\ne$PS (task 1, separate plan).

**Smearing = block-slice the full perambulator (NO new solves).** The distillation modes in $V$ are ordered by
Wilson eigenvalue, so keeping the first $n$ modes is the smearing projector $P_n=V_n V_n^\dagger$
($V_n=V[:,:n]$). A propagator leg from a source-vertex smeared at $m$ to a sink-vertex smeared at $n$ is the
$n\times m$ block of the full perambulator:
$$
\tau_{nm}(t',t)=V_n^\dagger(t')\,D_{ov}^{-1}(t',t)\,V_m(t)=\tau(t',t)[:n,:m].
$$
Position-propagator block (point machinery), smeared $n$ at sink $t'$, $m$ at source $t$:
$$
A_{nm}(t',t)=V(t')[:,:n]\;\tau(t',t)[:n,:m]\;V(t)[:,:m]^\dagger .
$$
Equal-time GW contact is improved in MODE space, $\tilde\tau_n(s,s)=\tau(s,s)[:n,:n]-\tfrac12 I_n$, so the
position contact is $-\tfrac12 P_n$ (NOT $-\tfrac12 I$); at $n=24$ (full, $V V^\dagger=I$) this reduces to the
existing $-\tfrac12 I$. Multi-smearing variational basis: Morningstar-Peardon hep-lat/9901004; distillation
Peardon 0905.2160; GPOF/Hankel Aubin-Orginos 1010.0202; GEVP Blossier 0902.1265.

## Operator basis
$\{\sigma^2_{00}(Y00),\,O_{2m}(\text{antipodal}),\,O_{1m}(\text{coincident split})\}$ (the full point set, PS
flavor = forward-improved legs) $\times$ smearings $\{4,8,16,24\}$ = up to **12 operators**. Correlator matrix
$C_{(i,n),(j,m)}(dt)=\langle O_i^{(n)}(t')\,O_j^{(m)}(t)\rangle_c$ includes all cross-smearing pairs (each leg
sliced by the smearings of the two vertices it connects). GEVP metric handles redundancy/rank.

## Chunks

### Chunk 1 -- smearing-aware point-correlator matrix (validate vs Nv=24)
Generalize the `fs_gevp_point`/`fs_channels_v2` point machinery so `AblkS` takes a per-vertex smearing level:
`AblkS_nm(ta,tb) = V[ta][:,:n] @ tau[ta,tb][:n,:m] @ V[tb][:,:m].conj().T`, equal-time contact `-0.5 * P_n`.
Build the 12x12 (3 ops x 4 smearings) connected matrix, PS (forward-improved) legs, config-parallel cache.
- VALIDATION: the $n=24$ sub-block must reproduce the existing `fs_channels_v2` 3x3 (PS=FS) to machine precision.
- Files: NEW `sigma2_smear_gevp_v2_claude.py` (reuses `fs_gevp_point_claude.make_config`/`perm_contrib_folded`
  with a smearing-level argument threaded through; cache `sigma2_smear_cache_claude/`).

### Chunk 2 -- GEVP (+ Hankel+rebase) on the expanded basis
Run the fixed-$t_0$ GEVP and the canonical Hankel+rebase (Dt=[0,2,4], rebase@4, T0=3) on the 12-op matrix;
compare state0/state1 (and any newly-resolved state) to the 3-op Nv=24 result. Binsize 10 jackknife. Does the
smearing basis sharpen $m_1$ / resolve a 3rd state cleanly?
- Files: extend the Hankel driver (`fs_channels_v2_hankel_claude.py` pattern) to read the 12-op cache.

### Chunk 3 -- (LATER) flavor-cross matrix
$\{\sigma_{PS}^2,\sigma_{FS}^2,\sigma_{FS}\sigma_{PS}\}$ with mixed furnishing (FS$\ne$PS). Separate plan; derive
the mixed-furnishing per-loop rule first ($-\tilde\tau$ on FS vertices vs $+\tilde\tau$ on PS; what survives of
Eq 5.5 $\langle S_4\rangle+\langle\tilde S_4\rangle$). Note $\sigma_{FS}\sigma_{PS}$ is parity-odd in continuum.

## Open questions
- Smearing levels: plan uses $\{4,8,16,24\}$ (24 = existing anchor). OK, or strictly $\{4,8,16\}$?
- Start on _v1 Nf2 gsq1.0 L1, 400 cfg (confirmed).

## CHUNK 1-2 RESULT (2026-09-10) -- smearing expansion is MARGINAL (mostly negative)
- Chunk 1 DONE: `sigma2_smear_gevp_v2_claude.py`, $n{=}24$ sub-block reproduces `fs_channels_v2` to $10^{-19}$.
  Fixed a per-config symmetrization bug ($O_{1m}$ time-split -> raw matrix asymmetric -> symmetrize at ENSEMBLE
  level). 12-op cache built (400 cfg, 12 threads ~24 min): `sigma2_smear_cache_claude/sigma2_smear_*_sm4-8-16-24_d1`.
- **Full 12-op is INDEFINITE** (metric eig $-0.6$ at $t{=}3$): mixing the 3 DIFFERENT point ops at TRUNCATED
  $N_v$ washes out the antipodal/coincident spatial structure -> the ops become redundant. Only per-operator
  smearing towers ($\sigma^2_{00}\times\{4,8,16,24\}$ etc.) are well-conditioned (PD). 3-ops@sm8 already indefinite.
- **Hankel on the smearing basis = over-parametrized garbage** (smearing is already a variational expansion;
  Hankel double-expands -> rank-deficient NaN). Use PLAIN pruned fixed-$t_0$ GEVP (`sigma2_smear_plaingevp_claude.py`).
- **Plain pruned GEVP**: descends from above (no fast plateau). Smearing augmentation (3-op@24 + $\sigma^2_{00}@\{8,16\}$)
  pulls state0 down marginally (t=5: 0.835 vs baseline 0.916) but gives NO new state and NO sharper two-meson
  ($m_1\sim0.62$ noisy either way). Physically: the two-meson needs the HIGH-mode antipodal $O_{2m}$ structure
  that smearing removes -> smearing helps the smooth ground, not the two-meson.
- **VERDICT: the 3-op Hankel+rebase (time-shift expansion) remains the best tool**; $N_v$-smearing expansion is
  marginal for this spectrum. Fig `figs/sigma2_smear_plaingevp_*`.
