# Adding the single-meson (m_PS) operator to the P+ GEVP basis -- impl plan

## Goal / physics
At the L2 **truncated** distillation basis (Nv=24 of 84), the single meson m_PS leaks into the sigma^2 (P+)
correlators and dominates the GEVP ground (verified: free complete Nv=84 -> two-meson 0.79; truncated Nv=24 ->
m_PS 0.40). The mechanism is NOT sigma3-hermiticity (the system lacks it); it is a truncation/smearing effect that
vanishes only at the complete basis. Rather than wait for Nv=84, **add an explicit single-meson operator sigma_00
(the ell=0 zero-momentum scalar bilinear) to the GEVP basis**, so the GEVP assigns the leaked m_PS to its own
state (via sigma_00) and the two-meson to the sigma^2 operators. NM: "I'll consider the sigma^2-m_PS channel" --
i.e. study the sigma_00 <-> sigma^2 mixing directly.

## Operators
- Single meson: sigma_00(t) = sum_x A_x Y00 psibar(x)psi(x) = ONE bilinear (2 fermions).
- sigma^2 ops (existing P+): {PP,FF} x {s2 (sigma^2_00), O2m (antipodal), O1m (coincident split)} = 4 fermions.
- Augmented basis (start minimal, extensible): {sigma_00, sigma^2_00(PP s2)} 2x2 -> the sigma^2-m_PS channel.
  Extend to {sigma_00} + the 6-op sigma^2 block (7x7) once the 2x2 is validated.

## Correlator matrix (all from the SAME L2 Nv=24 peram, same normalization)
- C_11 = <sigma_00(t) sigma_00(0)> : single-meson 2pt = single fermion loop through 2 vertices.
- C_22 = <sigma^2_00(t) sigma^2_00(0)> : the sigma^2 diagonal (the flavfac PP s2 four-point; reuse fs_gevp_point
  4-vertex machinery).
- C_12 = <sigma_00(t) sigma^2_00(0)> : the TRIANGLE (sink 1 bilinear, source 2 bilinears) -- single fermion loop
  through 3 vertices. C_21 = time-reverse.
All contractions use the plain forward AblkS = U tau U^dag (contact -1/2 on equal-time), tau ONLY (never tau_gw;
FS==PS collapse). NO sigma3 assumptions anywhere.

## Contractions (bespoke, since perm counts differ from the 4-vertex sigma^2)
- 2-vertex (single meson): C_11(dt) = -Tr[Phi(t) A(t,0) Phi(0) A(0,t)]  (one loop, (-1)^1 sign; A=AblkS spin/site
  blocks; Phi = area-Y00 mode vertex).
- 3-vertex (triangle): one fermion loop through {sink i @ t, src k,l @ 0}; sum the 2 orderings (i->k->l->i,
  i->l->k->i); flavfac for PP (all-PS -> factor per loop). Cross-check: reduces to the {1,1,1,1} ground_coupling
  form and to the flavfac cache C[0,0] for C_22.

## Files
- `sigma2_mPS_gevp_claude.py` (new): builds C_11, C_12, C_22 from the peram, assembles the augmented matrix,
  rebased GEVP, plots states with m_PS (0.3527 L2 / 0.3209 L1) and 2m_PS lines. Env: ENS, LREF, NVDIR, OPS, T0,
  REBT, NKEEP, BINSIZE. Reuses fs_gevp_point (AblkS, perm_contrib_folded, op_vspec).

## Open questions
1. Basis size: start 2x2 {sigma_00, s2}; extend to 7-op {sigma_00 + 6 sigma^2}? (default: build 2x2 first, show
   the channel, then extend.)
2. m_PS is 2-fermion, sigma^2 is 4-fermion -> the GEVP mixes operators of different fermion number; that is fine
   (they share the same 0++ states) but the C_12 triangle is the key new object. Validate C_12 != 0 and its sign.
