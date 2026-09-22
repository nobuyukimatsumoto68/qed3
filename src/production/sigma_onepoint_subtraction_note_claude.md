# Scalar one-point functions from data, and why the measured subtraction equals the contact

**Ensemble:** Nf2 $g^2{=}0.5$ L1, $a_t{=}0.2$, massless, 400 cfg. Distillation $N_v{=}24$ (exact), window $t_\text{win}{=}32$.
**Codes:** `one_point_claude.py` (one-point measurement), `two_meson_op1sub_claude.py` (contact vs measured-one-point subtraction), `diag_effmass_claude.py:20` (`diags_pair`).

## Building blocks (PS legs, $\tau=V^\dagger D_\text{ov}^{-1}V$)

Per window-time $s$, with a scalar equal-time subtraction $c$:
$$
D_S(s) = \mathrm{Tr}\big[\Phi(s)\,(\tau(s,s)-c\,I)\big],\qquad
D'_S(s) = \mathrm{Tr}\big[\Phi(s)(\tau(s,s)-c\,I)\,\Phi(s)(\tau(s,s)-c\,I)\big],
$$
$$
O(s) = 2\big(D_S(s)^2 + D'_S(s)\big)\quad\text{(composite }\sigma^2\text{ one-point)},\qquad \Phi(s)=V^\dagger(s)\,W_0\,V(s),\ \ W_0=A_x Y_{00}.
$$

## Measured one-point functions (config jackknife)

| | $\langle D_S\rangle$ | $\langle D'_S\rangle$ | $\langle O\rangle=\langle\sigma^2\rangle$ |
|---|---|---|---|
| raw ($c{=}0$) | $3.544908$ (std $1.7{\times}10^{-7}$) | $0.513751$ | $26.160$ |
| contact-subtracted ($c{=}\tfrac12$) | $0.000000$ (std $1.7{\times}10^{-7}$) | $-0.009848(3)$ | $-0.019700$ |

Two clean facts:

1. **$\langle D_S\rangle$ is *purely* the GW contact.** $D_S$ is config-independent to $10^{-7}$, and equals $\tfrac12\,\mathrm{Tr}\,\Phi = 3.5449$ (so $\mathrm{Tr}\,\Phi = 7.09$). After subtracting $c{=}\tfrac12$ it is $0.000000$ exactly: **no condensate, no non-contact correction** to the single-$\sigma$ tadpole. Consistent with no SSB.

2. **$\langle D'_S\rangle$ is NOT killed by the contact.** It leaves $-0.009848(3)$. The composite one-point $\langle O\rangle = 2\langle D_S^2+D'_S\rangle = -0.0197$ is *entirely* this non-contact $D'_S$ piece (since $D_S{=}0$ after the contact). This is the "additional correction" beyond the analytic $\tfrac12$.

## Subtracting the full measured one-point $\equiv$ the contact (for the effmass)

Replace the scalar $\tfrac12 I$ by the full data-measured equal-time one-point **matrix**
$$
M(s) = \langle\tau(s,s)\rangle_\text{cfg}\quad(N_v\times N_v),\qquad \tau(s,s)\to\tau(s,s)-M(s)\ \text{in every diagram.}
$$
Result (`two_meson_op1sub_claude.py`): the summed two-meson effmass is **unchanged to $<10^{-3}$** (e.g. $a_t m_\text{eff}(dt{=}17)=0.5177$ contact vs $0.5174$ matrix).

Reason: the contact already saturates the *linear* one-point. Numerically $\langle D_S\rangle_{c=1/2}=\mathrm{Tr}[\Phi(M-\tfrac12 I)]=0$, i.e. the measured one-point projected on the vertex is exactly $\tfrac12$:
$$
\mathrm{Tr}\big[\Phi\,(\langle\tau(s,s)\rangle-\tfrac12 I)\big]=0 .
$$
So every diagram carrying $\tau(s,s)$ *linearly* (the tadpoles $D_S$ in G, C, D, H, I, J, and the equal-time insertions inside the connected loops A, $V_S$, $S_S$) sees the same subtraction from $M$ as from $\tfrac12 I$. The residual $\langle D'_S\rangle=-0.00985$ is a **quadratic** (variance) effect — the normal-ordering of $\sigma^2$ — and does not move the linear one-point.

## Consequence for the sub-0++ plateau

The tadpole diagrams are already $\sim0$ after the contact ($|C_\text{conn}(dt{=}4)|\sim10^{-9}$; `diag_effmass_claude.py` table), so they do **not** cause the summed effmass ($\approx0.52$) to sit below the 0++ ($a_t m_{F^2}=0.616$) and below the clean two-$\sigma$ (E $=C_S^2\approx0.72$, B $=T_S\approx0.69$). The drag comes from the genuinely **connected** diagram $A=-S_S$ (single fermion loop visiting $s,s,t,t$), whose effmass falls through $0.78\to0.62\to0.53$ and carries intrinsic **single-$\sigma$** content by pinching. One-point subtraction cannot remove this — it is not a vacuum/tadpole artifact.

**Bottom line.** The one-point subtraction is now data-driven and complete (contact for the linear piece; $-0.00985$ measured for $D'_S$), but it does **not** raise the summed effmass to the two-$\sigma$ value. Isolating the true two-scalar state ($a_t m\approx0.70$, i.e. $\Delta/\Delta_A\approx2$, degenerate with the 0++) requires **diagonalizing** the operator basis (GEVP) or restricting to the manifestly two-$\sigma$ diagrams E, B — not a bound state below threshold.
