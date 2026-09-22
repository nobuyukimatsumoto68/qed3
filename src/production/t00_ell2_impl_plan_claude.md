# Stress tensor dimension: $\ell=2$ projection + singlet (disconnected) -- impl plan

## Physics / goal (NM 2026-09-21)
- $T_{\mu\nu}$ is a spin-2 primary: on $S^2\times\mathbb R$ its state is the $\ell=2$ multiplet at $E=3/R$.
  The $\ell=0,1$ components of $T_{00}$ are the charges ($H$, $K\pm P$) and create nothing; pair modes
  $(n_1,n_2)$ reach $\ell=0$ only for $n_1=n_2$ ($2/R,4/R,\dots$). So the old $\ell=0$ $O_H$ plateau
  (free L1 0.563, interacting L1 0.459) is the $(2,2)$ meson = {1,1,1,1}, NOT $T$; "Delta: 2.74 -> 3.16" is retracted
  as a statement about $T$.
- Operator: $O_H^{(m)} = \eta^\dagger W^{(m)}\xi + {\rm h.c.}$, $W^{(m)}_{ij}=W_{ij}\,Y_{2m}(\hat n_{\rm mid}(ij))$, real
  $Y_{2m}$ (keeps $W^{(m)}$ anti-hermitian, so $\Phi_{W^\dagger}=-\Phi_W$ and both halves of the connected 2pt are equal).
- Connected-only = flavor-ADJOINT spin-2 bilinear, NOT conserved: $\Delta_{\rm adj}=3+\gamma_{\rm adj}$, $\gamma=O(1/N_f)$.
- Singlet = $N_f C_{\rm conn}+N_f^2 C_{\rm disc}$: the disconnected piece is what restores the protected $\Delta_T=3$
  (gluonic intermediate states; $T^G$ itself only an optional extra GEVP operator).
- Per-config loop (correct equal-time legs: $\tau$ for the $W$ half, $1-\tau$ for the $W^\dagger$ half):
  $L_m(t) = -2\,{\rm Re}\,{\rm Tr}[\Phi^{(m)}(t)\,(\tau(t,t)-\tfrac12)]$; $C_{\rm disc}(dt)=\langle L_m(t+dt)L_m(t)\rangle-\langle L\rangle^2$.

## Files
- NEW `t00_ell2_claude.py` (driver: ELLS list, connected m-avg correlator, loops cached per config).
- reuse (read-only): `t00_stress_ham_interacting_claude.py` (build_W_gauge, read_config_sp, load_link_table),
  `geom_hopping_claude.py`, `distill_contract_claude.py`.
- update at end: `t00_hamiltonian_derivation_claude.md`, memory `project_t00_stress_tensor.md`.

## Chunks
- **A. Free $\ell=2$ (L1 Nv24, L2 Nv84; both complete).** Files: `t00_ell2_claude.py`.
  Expect $E_{1/2}+E_{3/2}$ ($\to 3/R=0.567$ in the continuum); reference lines from $\ell=0$: $\sigma$ 2pt $=2E_{1/2}$,
  $O_H^{\ell=0}=2E_{3/2}$. Also settles the $\ell=0$ misidentification (free L2 $\ell=0$: $0.567$ vs $2E_{3/2}(L2)$).
- **B. Interacting connected $\ell=2$ (adjoint), Nf2 g1.0 L1 then L2.** Files: same driver. Gives $\Delta_{\rm adj}$.
- **C. Singlet: add $C_{\rm disc}$ from the cached loops** (vacuum-subtracted, jackknife). Files: same driver
  (+ small analysis script if needed). Expect $\Delta_T=3$; noise-limited (loop variance const in $t$).
- **D. (optional) GEVP with $T^G$ / $O_T^{\ell=2}$.**

## Open questions
- L2 interacting is Nv24-of-84 TRUNCATED: $Y_{2m}$-weighted hop may be poorly represented; judge after chunk A/B.
- Corrections owed to the other threads: "<T_00 sigma^2>=0 exact" holds for the CONNECTED part only (tadpole leg fix).
