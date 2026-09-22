# The (2,2) single-meson operator via the ell=3/2 shell projector -- impl plan

## Physics / goal
Split the (2,2) single-meson excitation (Delta=4) from the two-meson (also Delta=4) inside the P+ sigma^2 GEVP
Delta=4 ~0.73 level. The (2,2) interpolator (from the {1,1,1,1} agent, verified free L1 -> 2E_{3/2}=0.556 clean):
Phi_(2,2)(t) = P_shell(t) = orthogonal projector onto the ell=3/2 eigenspace of the stored anti-hermitian
M(t) = tau(t,t) - 1/2 I_Nv (mode-space, equal-time contact-subtracted peram). NO new solve, NO Xi -- built from the
stored /peram/tau. Meson correlator C(t,s) = -Tr[P_shell(t) tau(t,s) P_shell(s) tau(s,t)] (all mode-space, Nv x Nv).

## Shell selection ({1,1,1,1}, verified)
- PER-TIMESLICE projector (unique orthogonal projector -> no eigvec phase/order ambiguity; tracks per-t fluctuation).
- RANK-based on |eig(M)|: |eig| ordering ell1/2=0.131(x4, ranks 1-4), ell3/2=0.078(x8, ranks 5-12),
  ell5/2=0.023(x12, ranks 13-24). Pick ranks 5-12 (0-indexed 4..11) = the 8-mode (2,2) shell. Count-stable interacting.
- Interacting caveat: degeneracy lifts (no rotational symmetry); if the (2,2) GEVP state looks contaminated, widen to
  ranks 1-12 (LO=0 HI=12, both low shells -> GEVP+sigma^2 split them) or use Phi22 as one variational op.

## snippet (author {1,1,1,1})
    def shell_projector(M, lo=4, hi=12):
        w, U = np.linalg.eigh(1j * M)          # iM hermitian; w real, U orthonormal
        sel = np.argsort(-np.abs(w))[lo:hi]     # ranks 5..12 = ell3/2 shell
        Us = U[:, sel]
        return Us @ Us.conj().T                 # Nv x Nv orthogonal projector

## Chunks
### Chunk 1 -- (2,2) DIAGONAL effmass (direct view of the state).  File: sigma2_shell22_claude.py (new)
Load /peram/tau (all windows), per-t Phi22 = shell_projector(tau[a,a]-1/2 I), C(dt)=translation-avg
-Tr[Phi22[t] tau[t,s] Phi22[s] tau[s,t]].  Jackknife effmass.  Interacting L2 -> expect a Delta=4 state (free 0.556;
the (2,2)).  Plot with m_PS, 2m_PS lines.  MODE_CONTACT irrelevant to the diagonal (contact is inside Phi22 via M).

### Chunk 2 -- add Phi22 as the 7th op to the P+ GEVP.  Files: sigma2_shell22_gevp_claude.py (new; models sigma2_mPS_gevp)
Basis {PP-s2,PP-O2m,PP-O1m,FF-*,Phi22} (or start minimal {sigma^2_s2, Phi22}).  Need the CROSS
<Phi22(t) sigma^2(0)> (triangle: Phi22 mode-space sink vertex + 2 sigma^2 vertices) -- reuse perm_contrib_n with a
mode-space Phi22 sink.  GEVP -> the (2,2) pulls to its own state, two-meson stays in the sigma^2 block.

## Open
- Is Phi22 sink cross with the O2m/O1m GEOMETRY (antipodal / split) well-defined with a mode-space projector vertex?
  s2 (local) is straightforward; O2m/O1m carry position-space geometry -- may need the position-space triangle.
  Start with {sigma^2_s2, Phi22} (both expressible), then extend.
