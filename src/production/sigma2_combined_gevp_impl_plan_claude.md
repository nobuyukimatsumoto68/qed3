# Combined (2,2)-shell + sigma^2 two-meson GEVP -- impl plan

## Goal (NM)
"Three states around 2m_PS" -- the (2,2) single meson (~0.65), the two-meson (~0.72), an excited (~0.80). Build ONE
GEVP over the (2,2) shell single-meson operator(s) AND all the sigma^2 two-meson operators, keep NKEEP=4, so all of
them are resolved and orthogonalized together.

## Operators
- SHELL single-meson (PS): P_shell(t) on M-eigenspace |eig| rank windows.  Default 3: ell1/2=[0,4], ell3/2=[4,12]
  (the (2,2)), ell5/2=[12,24].  Vertex nonlocal K_shell(t) = U(t) P_shell(t) U(t)^H (2Ns x 2Ns).  (ell1/2 = m_PS
  clearer.)
- sigma^2 two-meson (PP; PP==FF so drop FF redundancy): s2 (op0), O2m (op1, antipode), O1m (op2).  op_vspec geometry.

## Correlator blocks (per config, translation-avg over source s, sink t=s+dt)
- shell x shell (mode space):  C = -Tr[ P_shell(t) tau(t,s) P_shell(s) tau(s,t) ].
- sigma^2 x sigma^2 (PP):  the flavfac-PP 4-vertex, perm_contrib_folded + FFPP (as flavorgeom / sigma2_mPS C_22).
- CROSS shell x sigma^2 (the new piece): 3-vertex triangle, ONE fermion loop through {shell sink @ t, src a @ s,
  src b @ s}.  = -Tr[ K_shell(t) A(t,s) Wa A(s,s)_sub Wb A(s,t) ] + (a<->b), A = AblkS (position blocks, MODE_CONTACT
  on the equal-time A(s,s)), Wa/Wb/antipode from op_vspec(op).  Implemented with the shell as a MATRIX vertex
  (2 site-letters + K_shell factor) in a bespoke einsum, or by folding into a blob B(s,s)=A(s,t)K_shell(t)A(t,s).
  VALIDATION: with K_shell = identity (P_shell=I_Nv, i.e. window [0,Nv]) the cross must reduce to the local
  <psibar psi | sigma^2> triangle == sigma2_mPS_gevp C_12 (for op0/s2).

## GEVP
Assemble (3 shell + 3 sigma^2) = 6x6 (extensible), no-Hankel (OFFSETS=0) rebased GEVP, NKEEP=4.  Expect:
state0 m_PS (from ell1/2), state1 (2,2)~0.65 (ell3/2), state2 two-meson~0.72 (sigma^2), state3 excited~0.80.

## Files
- sigma2_combined_gevp_claude.py (new).  Reuses G=fs_gevp_point (make_config_win, op_vspec, perm_contrib_folded),
  MG=sigma2_mPS_gevp (perm_contrib_n, FFPP), shell_projector ({1,1,1,1}).  MODE_CONTACT=1.

## Validate order
1. K_shell=I cross == sigma2_mPS C_12 (s2).  2. FREE L1: combined GEVP resolves m_PS 0.378, (2,2) 0.556, two-meson
0.756 as distinct states.  3. Interacting L2, NKEEP=4.
