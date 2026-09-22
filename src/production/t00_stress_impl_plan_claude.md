# Fermionic stress-tensor $T_{00}$ from the distillation perambulator -- impl plan

## Goal / physics
Measure the energy-density component of the fermionic stress tensor, $T_{00}$, on $S^2\times R$ using the
EXISTING perambulator (no new solves). $T_{00}$ is a fermion bilinear with a temporal derivative, so its
two-point is a single distillation loop -- the same machinery as the $\sigma=\bar\psi\psi$ two-point, only the
VERTEX changes. Motivation: identify the excited $0^{++}$ ("{1,1,1,1}") state seen in the two-meson correlators
-- $T_{00}$ is the universal $\Delta=3$, $\ell=0$, $0^{++}$ operator (stress tensor is protected, $\Delta=d=3$),
in the SAME channel as $\sigma$/$\sigma^2$, so it can mix with them. Decisive check is the free limit, where the
stress tensor sits at exactly $\Delta=3$.

## Operator (interpolator, 2-component $S^2\times R$)
Temporal gamma is $\sigma_3$ (the "s3/temporal" channel of the axial/vector trio; the free on-axis propagator
carries the $\sigma_3$ structure). The $\ell=0$ energy-density interpolator:
$$
\mathcal{O}_{T_{00}}(t) = \sum_x w_x\; \bar\psi(x,t)\,\sigma_3\,(\partial_0^{\rm sym}\psi)(x,t),
\qquad w_x = A_x\,Y_{00},\quad \partial_0^{\rm sym}\psi(x,t)=\tfrac12\big[\psi(x,t{+}1)-\psi(x,t{-}1)\big],
$$
with $A_x$ = `dual_areas` and $Y_{00}=1/\sqrt{4\pi}$ (the $\ell=0$ projection, exactly the $\sigma$ weight).
- Quantum numbers: $0^{++}$ (parity-even), $\ell=0$, $\Delta=3$ (free: $\Delta_\psi=1$ each $+1$ derivative) --
  same channel as $\sigma$ (Δ=2) and $\sigma^2$, so it mixes with the {1,1,1,1} level.
- This is an INTERPOLATOR for the state, NOT the exactly-conserved lattice stress tensor: the properly
  conserved overlap $T_{\mu\nu}$ needs GW point-splitting (like the conserved current's $(1-D_{ov})$), but for a
  STATE energy the naive $\bar\psi\sigma_3\partial_0\psi$ is enough (we want $E$, not the Ward identity).
- Traceless caveat: the improved $T_{00}$ subtracts $\tfrac1d\delta_{00}\,\bar\psi(\gamma\!\cdot\!D)\psi$, but
  $\gamma\!\cdot\!D\,\psi$ is the Dirac EOM ($\to 0$ massless), so the trace piece is EOM-null -- negligible for
  the interpolator; keep the naive form for chunk 1.
- Contamination caveat: $\bar\psi\,\partial_0\psi$ (no $\sigma_3$) is the $\partial_t$-descendant of $\sigma$ and
  rides the $\sigma$ state (E unchanged), so keep the $\sigma_3$ so the operator projects the genuine $\Delta=3$
  primary; the GEVP with $\sigma$ (chunk 3) cleanly separates any residual $\sigma$-tower leakage.

## Contraction (connected single loop; reuse AblkS)
$\psi(x,t')\to$ position-space propagator block $S(x,t';y,t)=$ `AblkS(t',t)[x,:,y,:]` (the improved block
$V(t')\tau(t',t)V(t)^\dagger$ from `fs_gevp_point_claude.make_config`). The connected two-point:
$$
C_{T_{00}}(t_0,t) = -\sum_{x,y} w_x w_y\;
\mathrm{Tr}_{\rm spin}\big[\,\sigma_3\, \mathbf{D}_0 S(x,t;y,t_0)\, \sigma_3\, \mathbf{D}_0 S(y,t_0;x,t)\,\big],
$$
where $\mathbf{D}_0 S$ denotes the symmetric $\pm1$ temporal finite difference applied on BOTH the source and
sink time arguments (four AblkS blocks at shifted times $t\pm1$, $t_0\pm1$, assembled as $\tfrac12(\cdots)$).
Translation-average over $t_0$ within the window (as the $\sigma$ two-point does). Fold $C(dt)=C(-dt)$.
- Vacuum/disconnected: subtract $\langle\mathcal{O}\rangle^2$ (the self-contracted single-loop piece); check
  whether $\langle\mathcal{O}_{T_{00}}\rangle=0$ by symmetry (expected, temporal-parity odd insertion).
- CONTACT: AblkS carries the $-\tfrac12 I$ improvement at equal time; the $\partial_0^{\rm sym}$ never uses the
  equal-time block ($t'\ne t\pm1$ vs $t$), so no contact ambiguity in $\mathbf{D}_0 S$. Confirm in code.

## Validation targets (free limit, decisive)
Free peram `data_free/distill_Nv24/peram.0.h5` (L1, Nv=24 COMPLETE = exact), single config, no jackknife.
- $m_\sigma = 2/R = 0.378$ (lattice $a_t m$) sets $1/R=0.189$.
- Expect a clean $T_{00}$ plateau at $E = 3/R = \mathbf{0.567}$ (lattice), DISTINCT from $\sigma$ (0.378) and
  from the two-meson ground $4/R=0.756$. Landing at 0.567 validates the vertex + the $\Delta=3$ ID.

## Chunks
### Chunk 1 -- free-limit $T_{00}$ two-point + validation
Build $C_{T_{00}}(dt)$ on the free peram, effmass, confirm plateau at 0.567. Also print $\langle\mathcal{O}\rangle$
(expect $\approx 0$).
Files: NEW `t00_stress_claude.py` (reuse `fs_gevp_point_claude.make_config` -> AblkS, `distill_contract_claude`
dual_areas + Y00; single-config/free guard as in the sigma2 scripts). Fig `figs/t00_stress_free_*_claude.png`.

### Chunk 2 -- interacting $T_{00}$ (Nf2 gsq1.0 L1)
Same driver on `distill_Nv24` (400 cfg) and/or `distill_Nv24_v2` (1000 cfg, nsrc2 via load_peram_windows),
jackknife bin10. Effmass -> $\Delta_{T_{00}}$; compare to $\sigma$ and to the {1,1,1,1} level.
Files: extend `t00_stress_claude.py` (ENS/NVDIR env).

### Chunk 3 -- coupled $\{T_{00}, \sigma^2, \ldots\}$ GEVP (the mixing)
Cross-correlate $T_{00}$ with $\sigma^2$ (and the flavor/geom two-meson ops); block-Hankel+rebase GEVP. Read
whether the {1,1,1,1} state carries stress-tensor character (nonzero $\langle T_{00}\,\sigma^2\rangle$ at that
level = the identification).
Files: NEW `t00_sigma2_gevp_claude.py` (reuse the hankel_reb path + the flavor cache).

## Gluonic piece / the CONSERVED stress tensor (hand-off decision)
Do we need the gauge part $T^{\rm gauge}_{\mu\nu}=F_{\mu\alpha}F_\nu{}^\alpha-\tfrac14\delta_{\mu\nu}F^2$?
- The physical conserved tensor is the SUM $T_{\mu\nu}=T^{\rm ferm}_{\mu\nu}+T^{\rm gauge}_{\mu\nu}$; only the
  sum is conserved and $\Delta=3$-PROTECTED. Separately there are TWO $\Delta=3$ spin-2 operators (ferm, gauge):
  one linear combination is the conserved $T$ ($\Delta=3$ exact), the ORTHOGONAL combination is non-conserved and
  picks up an anomalous dimension in the interacting theory. So the fermionic $T_{00}$ alone overlaps BOTH -- a
  single-operator effmass sees a mixture, not the clean protected $\Delta=3$.
- FREE limit (chunk 1): gauge and fermions DECOUPLE, so the fermionic $T_{00}$ alone lands at exactly $\Delta=3$
  (0.567). The gluonic piece is NOT needed to validate the vertex or to see the free $\Delta=3$ state.
- The immediate {1,1,1,1} question is a FERMIONIC-sector question: $\sigma^2$ is purely fermionic, so the
  fermionic $T_{00}$ is the DIRECT overlap probe; the gluonic $T_{00}$ couples to $\sigma^2$ only through the
  gauge interaction (subleading). So chunks 1-3 (fermionic only) answer "does {1,1,1,1} carry stress-tensor
  character" on their own.
- To identify the state as THE conserved stress tensor in the INTERACTING theory, YES -- add the gluonic piece
  and do a joint $\{\sigma^2,\,T^{\rm ferm}_{00},\,T^{\rm gauge}_{00}\}$ GEVP; the conserved combination is the
  one that stays $\Delta=3$ across couplings.
- RECOMMENDATION: build the fermionic $T_{00}$ first (chunks 1-3, self-contained, reuses the perambulator). The
  gluonic $T^{\rm gauge}_{00}$ = traceless energy density of the flow-smeared field strength is a GLUE-SECTOR
  observable (same $F_{\mu\nu}$ / Wilson-flow machinery as the $F^2$ glueball) -> HAND OFF to the glue agent
  (qed3-a7); merge into the joint GEVP as chunk 4 for the final conserved-$T$ identification.

## Open questions -- RESOLVED with NM (2026-09-16)
1. $\gamma_0 = \sigma_3 = \mathrm{diag}(1,-1)$. CONFIRMED (temporal Pauli, s3 channel / free on-axis structure).
2. Stencil: START symmetric $\tfrac12(\psi_{t+1}-\psi_{t-1})$; may add one-sided/improved LATER for comparison.
3. START naive (no traceless subtraction); may check improved/traceless LATER.
4. Basis: free single config (chunk 1) -> Nf2 400 cfg -> nsrc2 1000. OK.
5. Vacuum: print $\langle\mathcal{O}_{T_{00}}\rangle$ (expect $\approx 0$ temporal-parity); subtract
   $\langle\mathcal{O}\rangle^2$ only if it comes out nonzero.
6. CONTACT (new, raised by stress-tensor agent): at $dt=1,s_x={-}1$ the block is the equal-time
   `AblkS(s,s)`. NM: "always better to subtract the contact term IN the perambulator" -- and `AblkS` already
   does this (`A = A - 0.5*Iv` at `ta==tb`). RESOLUTION: use `AblkS` as-is (contact already subtracted in
   the perambulator block); no special-casing of the $dt=1$ term.

### Contraction as coded (window-relative times; sink at $s{+}dt$, source at $s$, translation-avg over $s$)
$$
C_{T_{00}}(dt) = -\tfrac14 \sum_{s_x,s_y=\pm1} s_x s_y \sum_{x,y} w_x w_y\,
\mathrm{Tr}_{\rm spin}\!\big[\sigma_3\,\mathrm{AblkS}(s{+}dt{+}s_x,\,s)[x,:,y,:]\;
\sigma_3\,\mathrm{AblkS}(s{+}s_y,\,s{+}dt)[y,:,x,:]\big]
$$
$\partial_0^{\rm sym}$ differences the FIRST (out) time index of each block. $w_x = A_x Y_{00}$. Fold $C(dt)=C(-dt)$.
Single-config free guard: BINSIZE=1, `nb<2 -> em_err=0` (as in the sigma2 scripts). Effmass = $\log[C(dt)/C(dt{+}1)]$.

Refs: distillation Peardon 0905.2160; derivative/displaced distillation operators (Peardon-Edwards spectroscopy);
stress tensor $\Delta=d=3$ universal; free $S^2\times R$ propagator (Eq C.28, cont_prop) for the shell check;
sigma^2 machinery `fs_gevp_point_claude.py`, `sigma2_flavorgeom_full_v2_claude.py`; mixing context CP 1603.05582.
