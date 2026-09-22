# Extend interacting VSH sp to ell=2 -- impl plan

Add the ell=2 VSH tower (electric $\Phi_2$, magnetic $\Psi_2$, tp radial ell=2) to the interacting sp
pipeline, on the SAME perams / same site-level $f^{ab}(n_1,n_2)$ -- only the angular weights change.

## Dimension predictions (Delta=2; tp&$\Phi_\ell$ dim $\Delta+\ell-1$, $\Psi_\ell$ dim $\Delta+\ell$)
- electric $\Phi_2$ dim $\Delta+1=3$  -> should match tp ell=2 (H) = 0.5075(16) and my $\Psi_1$.
- magnetic $\Psi_2$ dim $\Delta+2=4$  -> should match tp ell=3 = 0.6001(55).
- tp radial ell=2 (my cross-check) -> should match d4 tp ell=2 = 0.5075(16).
L2 (42 sites) is where ell=2 is meaningful; L1 (12 sites) ell=2 ALIASES hard (H is 5-dim on icosahedron) -- run it but flag unreliable.

## Real $Y_{2m}$ and pole-safe tangent gradients (convention: grad[m]=(df/dtheta, (1/sin t) df/dphi), same as ell=1)
Let $A=\tfrac14\sqrt{5/\pi}$, $B=\tfrac12\sqrt{15/\pi}$, $C=\tfrac14\sqrt{15/\pi}$; $s=\sin t$, $c=\cos t$.
- m=0:  Y = A(3c^2-1);            grad=(-6A s c, 0)
- m=1:  Y = B s c cos(ph);        grad=(B cos(2t) cos(ph), -B c sin(ph))
- m=-1: Y = B s c sin(ph);        grad=(B cos(2t) sin(ph),  B c cos(ph))
- m=2:  Y = C s^2 cos(2ph);       grad=(C sin(2t) cos(2ph), -2C s sin(2ph))
- m=-2: Y = C s^2 sin(2ph);       grad=(C sin(2t) sin(2ph),  2C s cos(2ph))
All phi-components finite at the poles (the 1/sin t cancels), same as ell=1. The overall $1/\sqrt{\ell(\ell+1)}$
VSH normalization is a constant per ell -> cancels in the effmass, so it is dropped (as in the ell=1 code).
Then electric WE={1:grad_theta,2:grad_phi}; magnetic WM={1:-grad_phi,2:grad_theta}; tp uses the scalar Y value.

## Files
- Edit `interacting_vsh_L2_claude.py` and `interacting_vsh_L1_claude.py` (mine): add `ylm2_weights`, an `ELL`
  env switch (default 1) selecting weights + m-list ((-1..1) or (-2..2)), output names carry `_ell{ELL}`.
- Run ELL=2 for L2 (meaningful) and L1 (flagged aliased). Compare $\Psi_2$ to tp ell=3, $\Phi_2$ to tp ell=2.
