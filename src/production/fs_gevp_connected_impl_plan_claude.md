# FS connected GEVP with FULL-connected matrix elements (unified brute-force Wick)

Goal: build the FS 2-op connected GEVP $\{\sigma_{FS}^2,\ O_{\sigma\sigma}^{FS}\}$ with EVERY matrix element
computed from the SAME full-connected contraction rule (A,B,C,D,E,G), fixing the earlier inconsistency where
the $O_{2\sigma}$ sector (its diagonal AND its cross with $\sigma^2$) was E-only ($2M\!\cdot\!M$) while
$\sigma^2\sigma^2$ was full A-G. NM (2026-09-09): do FS first, establish the method, then redo PS.

## Operators (REVISED 2026-09-09, NM: geometry-optimized two-meson / one-meson interpolators)

Keep $\sigma^2_{00}$ (Y00) as the PRIMARY target; the two auxiliary interpolators are built by placing the two
$\sigma$ densities at chosen spatial/temporal locations so the GEOMETRY projects onto the wanted state (no
artificial diagram dropping -- full-connected correlators, the placement does the projection):

- **Primary** $\sigma^2_{00}(\tau)=\big(\sum_x A_x Y_{00}\,\sigma(x,\tau)\big)^2$ -- Y00-projected, vertices both at $\tau$.
- **Two-meson-optimal** $O_{2m}(\tau)=\sum_x A_x\,\sigma(x,\tau)\,\sigma(P(x),\tau)$ -- the two densities at
  ANTIPODAL points $x,P(x)$, EQUAL time. The one-meson contamination needs an equal-time meson to cross
  $x\to P(x)$, suppressed by $e^{-m_{PS}\,d(x,P(x))}$ with $d=\pi$ (maximal) $\Rightarrow$ one-meson killed.
- **One-meson-optimal** $O_{1m}(\tau)=\sum_x A_x\,\sigma(x,\tau)\,\sigma(x,\tau{+}\delta)$ -- SAME spatial point,
  time-split $\delta$ $\Rightarrow$ meson propagates a short time, no spatial suppression $\Rightarrow$ one-meson enhanced.

$P(x)$ = antipodal map, EXACT on L1 (12 sites, involution, $d=\pi$, err $10^{-16}$): P=[6,11,3,2,5,4,0,8,7,10,9,1].
$A_x$ = dual area. Point vertex $\Phi_x(\tau)=V^\dagger(\tau)\,P_x\,V(\tau)$, $P_x$ = projector onto site $x$'s
NS=2 spinor dofs (index $2x,2x{+}1$). Note $\Phi_{00}=Y_{00}\sum_x A_x\Phi_x$ (the Y00 vertex is the area-sum of point vertices).

### Efficient computation (avoid the naive $x,y$ double sum)
Form the position-projected propagator $Q(t,s)=V(t)\,\tau(t,s)\,V^\dagger(s)$ (a $2N_s\times2N_s$ matrix; $24{\times}24$
at L1). Then the point meson propagator $M_{xy}(t,s)=-\mathrm{Tr}_{\rm spin}[\,Q(t,s)_{xy}\,Q(s,t)_{yx}\,]$ (2x2
spin blocks), point tadpole $D_S(x,t)=\mathrm{Tr}_{\rm spin}[(\text{tt})(t,t)_{xx}]$, etc. The antipodal/coincident
sums become block sums over sites -- cheap. Build all connected diagrams (A-G) from these blocks (S-part leg
$\tau$, $\tilde S$-part leg $-\tau'$; FS = sum). This replaces the per-perm brute force for the point operators.

## GEVP basis
$\{\sigma^2_{00},\ O_{2m},\ O_{1m}\}$ (FS), full-connected, no identity. $O_{2m}$ overlaps the two-meson cleanly,
$O_{1m}$ the one-meson; $\sigma^2_{00}$ is the target whose content the GEVP resolves.

## (superseded) earlier 2-op set
- $\sigma^2(\tau)$: two $\sigma_{FS}$ vertices, both at time $\tau$.
- $O_{\sigma\sigma}(\tau)=\sigma_{FS,00}(\tau)\sigma_{FS,00}(\tau{+}\delta)$: Y00 both, time-split -- the earlier
  (Y00) two-$\sigma$; kept for reference. The point-based $O_{2m},O_{1m}$ above replace it.

## Unified connected matrix element (brute-force Wick)
For $\langle O_a(t)\,O_b(s)\rangle$, list the 4 vertex times = sink$(a,t)$ ++ source$(b,s)$, with
sink$(\sigma^2,t)=[t,t]$, sink$(O_{\sigma\sigma},t)=[t,t{+}\delta]$, source likewise at $s$. Then
$$
C_{ab}(t,s)=\sum_{\pi\in S_4}(-1)^{\#\text{cyc}(\pi)}\!\!\prod_{\text{cyc}}\mathrm{Tr}\big[\Phi\,\text{leg}\,\Phi\,\text{leg}\cdots\big],
$$
KEEP only permutations whose cycles BRIDGE sink$\{0,1\}\leftrightarrow$source$\{2,3\}$ (= connected A,B,C,D,E,G;
drop disconnected F,H,I,J). Contact on coincident legs: $S$-part $\tau(a,a)-\tfrac12 I$; $\tilde S$-part
$-\tfrac12(\tau'+\tau)(a,a)$. FS $=$ $S$-part(leg $\tau$) $+$ $\tilde S$-part(leg $-\tau'$). Same rule for
ALL of $C_{00},C_{01},C_{10},C_{11}$ $\Rightarrow$ consistent GEVP. Off-diagonal symmetrized $C_{01}\!\leftarrow\!\tfrac12(C_{01}{+}C_{10})$.

Matrix-element vertex lists (t = sink, s = source, $t=s{+}dt$):
- $C_{00}=\langle\sigma^2\sigma^2\rangle$: $[t,t,s,s]$
- $C_{01}=\langle\sigma^2 O_{\sigma\sigma}\rangle$: $[t,t,s,s{+}\delta]$
- $C_{10}=\langle O_{\sigma\sigma}\sigma^2\rangle$: $[t,t{+}\delta,s,s]$
- $C_{11}=\langle O_{\sigma\sigma}O_{\sigma\sigma}\rangle$: $[t,t{+}\delta,s,s{+}\delta]$

## Chunks
1. **Machinery + cross-check** (`fs_gevp_connected_claude.py`): module-level `wick_connected(times, Phi, legfn)`;
   verify $C_{00}$ (brute, $S$-part only, $\times2$) reproduces the PS `diag_corr_linear` connected A-G sum, and
   that $C_{01}$'s two-meson part matches the old $2M\!\cdot\!M$ up to the extra A,B,C,D,G. Report.
2. **FS 2x2 connected GEVP**: assemble $C(dt)$ (config-avg, binsize-10 jackknife), plain GEVP (metric $t_0$),
   effmass for both states; NO identity (vacuum-free by connected selection). Linear/log + effmass plot.
3. **Hankel+rebase** on the FS matrix (reuse `hankel_rebase_scan` machinery) if the plain GEVP is contaminated.
4. **THEN PS**: rerun `gevp_twosigma_interacting` sector with the full-connected $O_{2\sigma}$ (replace the
   E-only `Css`/`C2s`), re-extract $m_0,m_1$; compare to the old $0.46/0.6145$ (the sub-threshold claim).

## Cost
Brute-force = 24 perms/element. 4 elements x 2 legs x ~30 sources x ~16 dt x 400 cfg. ~0.8 s/config/element
measured -> ~15-20 min full. Test on a few configs first; run full detached if needed. OMP=4.

Refs: distillation Peardon 0905.2160; GEVP Blossier 0902.1265; mixing Chester-Pufu 1603.05582.
