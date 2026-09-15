# C++ data-dumper for the FS point-operator connected GEVP matrix (spec)

Fast replacement for the Python `fs_gevp_point_claude.py` heavy compute, for the production sweep
(gsq0.5/1.0/1.5 x Nf2/4/6 x L1/L2, PS and FS). The Python driver is the VALIDATED reference (chunks 1-2
exact to machine precision); the C++ must reproduce its per-config matrix $C_{ab}(t)$ bit-for-bit (up to
FP), then Python only reads the h5 and does the jackknife GEVP.

> STATUS: write the CODE only AFTER the Python FS L1 gsq1.0 result confirms the 3-op method resolves the
> spectrum (operators sufficiently independent; metric well-conditioned at full stats). Otherwise the
> operator set may change and the C++ would be labor in vain (NM's principle). Spec is safe to fix now.

## Physics / operators (identical to the Python driver + `fs_gevp_connected_impl_plan_claude.md`)
- $\sigma^2_{00}$ (Y00, separable kernel $w_iw_j$, $w=A_x Y_{00}$), $O_{2m}$ (antipodal $A_x,\ x{\to}P(x)$,
  equal time), $O_{1m}$ (coincident $A_x$, time-split $\delta$). All FS; full-connected (keep A-G, drop F-J
  = the bridging permutations); FS $=$ $S$-part (leg $\tau$) $+$ $\tilde S$-part (leg $-\tau'$).
- Contact: $S$ equal-time diagonal $\tau(t,t)-\tfrac12 I$; $\tilde S$ equal-time $-\tfrac12(\tau'+\tau)(t,t)$.

## Method (mirror the validated Python)
1. Load per config from `data_<ENS>/<NVDIR>/peram.<k>.h5`: `V` (twin,Nv,2Ns), `peram/tau`, `peram/tau_gw`
   (twin,twin,Nv,Nv) complex, meta tsrc0/twin. Use HighFive + Eigen (as in `glue_gevp_analysis_claude`).
2. Position propagator $A_{\rm leg}(t,s) = U_t\,{\rm leg}(t,s)\,U_s^\dagger$, $U_t=V[t]^{\rm T}$ (2Ns x Nv);
   store as Eigen `MatrixXcd` (2Ns x 2Ns). Site 2x2 spin blocks $A(t,s)_{xy}$.
3. Brute-force connected sum over the 20 bridging permutations of the 4 vertices (sink v0,v1; source v2,v3),
   value $=(-1)^{\#cyc}\prod_{\rm cyc}{\rm Tr}_{\rm spin}[\ldots]$, with kernel FOLDING (no dense site sums):
   - $\sigma^2_{00}$: weight each vertex site by $w$ (independent).
   - $O_{1m}$: the two vertices share a site (sum $\sum_x$), weight $A_x$.
   - $O_{2m}$: v1 site $=P(v0)$ (antipode-reindex the incident blocks), weight $A_x$.
   In C++ the site sums are plain loops (fast) -- no einsum-path issue; loop order: perm -> cycle -> sites.
4. Translation-average over source $s_0$ in the window; produce $C_{ab}(dt)$, $a,b\in\{0,1,2\}$, $dt=0..$twin-1,
   both legs summed. Symmetrize $C_{ab}=C_{ba}$.
5. OpenMP parallel over configs (or over $dt$); each config independent.

## Output
`data_<ENS>/fs_gevp_point_<NVDIR>/corr.<k>.h5` with dataset `C` shape (3,3,twin) real (+ meta SPLIT, tsrc0,
twin, op labels). Python analysis (`fs_gevp_point_claude.py` in an h5-read mode) globs these, does binsize
jackknife GEVP + Hankel. One file per config -> resumable, parallel-safe.

## Validation gate (mandatory before trusting the sweep)
On one config, C++ `C[a,b,dt]` must equal the Python `matrix_one_config` output to ~1e-10 (and the Python
Y00 element already equals `diags_pair` connected). Print the max abs diff.

## Build / handoff
`fs_gevp_point_dump_claude.cu` (or .cpp; no CUDA needed -- CPU Eigen). Compile with the project's HighFive/
Eigen includes (see `glue_gevp_analysis_claude.o` build guard). Handoff `tmp_claude.sh`: build-if-stale,
`OMP_NUM_THREADS=<cores>`, loop configs, tee to `*_claude.log`. NO rm/kill in the script.

## Files
- `fs_gevp_point_dump_claude.cu` (the dumper), handoff script, this plan.
- Python reference (validated): `fs_gevp_point_claude.py`; geometry/antipode from `distill_contract_claude.py`.
Refs: distillation Peardon 0905.2160; GEVP Blossier 0902.1265; mixing Chester-Pufu 1603.05582.
