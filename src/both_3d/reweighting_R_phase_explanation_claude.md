# Why the reweighting factor $R$ is a pure phase

Companion note to `reweighting_R_claude.cu` / `reweighting_R_impl_plan_claude.md`.
Reference: `qed3int_v2-14.pdf`, Eqs. (2.4)-(2.5).

## The factor

For the parity-breaking case the mass is purely imaginary, $m_P = i\mu$ with $\mu\in\mathbb{R}$.
The code computes (Eq. 2.5)

$$
R \;=\; \left(\prod_i \frac{(1-m_P)\lambda_i + m_P}{\lambda_i + m_P}\right)^{\!*},
$$

where $\{\lambda_i\}$ are the eigenvalues of the **massless** overlap $D_\text{ov}$. The conjugation
does not affect $|R|$, so "$R$ is a phase factor" means

$$
|R|^2 \;=\; \prod_i \frac{\big|(1-m_P)\lambda_i + m_P\big|^2}{\big|\lambda_i + m_P\big|^2} \;=\; 1 .
$$

## The key fact: the Ginsparg-Wilson circle

The overlap in this code is normalized as $D_\text{ov} = 1 + V$ with $V$ unitary (see
`overlap_wmass_claude.h:371`: `d_res = C*(D_W/lambda_max)*rational(H) v + 1*v`, the rational part
approximating $\gamma\,\mathrm{sgn}(H_W)$). Hence every eigenvalue of $D_\text{ov}$ lies on the
**GW circle of radius 1 centered at 1**,

$$
|\lambda_i - 1| \;=\; 1 \qquad\Longleftrightarrow\qquad \lambda_i + \bar\lambda_i = |\lambda_i|^2 .
$$

This is the eigenvalue form of the GW relation $D_\text{ov} + D_\text{ov}^\dagger = D_\text{ov}^\dagger D_\text{ov}$
and holds on **any** gauge background, not just $U=1$.

## The cancellation (per eigenvalue, exact)

Write $m_P = i\mu$ and compute numerator and denominator magnitudes for a single eigenvalue $\lambda$.

Denominator:

$$
|\lambda + i\mu|^2
 = (\lambda + i\mu)(\bar\lambda - i\mu)
 = |\lambda|^2 + i\mu(\bar\lambda - \lambda) + \mu^2
 = |\lambda|^2 + 2\mu\,\mathrm{Im}\,\lambda + \mu^2 .
$$

Numerator: let $w = (1-i\mu)\lambda + i\mu = \lambda + i\mu(1-\lambda)$. Then

$$
|w|^2
 = |\lambda|^2 + i\mu\big[\bar\lambda(1-\lambda) - \lambda(1-\bar\lambda)\big] + \mu^2|1-\lambda|^2
 = |\lambda|^2 + 2\mu\,\mathrm{Im}\,\lambda + \mu^2\,|1-\lambda|^2 .
$$

(The cross term reduces to $\bar\lambda-\lambda = -2i\,\mathrm{Im}\,\lambda$, identical to the denominator's.)

The numerator and denominator differ **only** in the last term: $\mu^2|1-\lambda|^2$ vs. $\mu^2$.
On the GW circle $|1-\lambda| = 1$, so

$$
\boxed{\;|w|^2 = |\lambda + i\mu|^2 \quad\text{for every } \lambda \text{ with } |\lambda-1|=1\;}
$$

Therefore the ratio has modulus 1 **eigenvalue by eigenvalue**, and $|R| = 1$ exactly. No pairing of
$\lambda$ with $\bar\lambda$ is even needed; the GW circle alone does it. (Conjugation symmetry of the
spectrum, from $\gamma_5$-hermiticity $D_\text{ov}^\dagger=\gamma_5 D_\text{ov}\gamma_5$, additionally
makes $R$ real $\times$ phase consistent, but is not required for $|R|=1$.)

Two remarks:
- The result is **special to purely imaginary $m_P$** (the parity case). For a general complex $m_P$ the
  two extra terms no longer match and $|R|\neq 1$ in general. A purely real $m_F$ would also not give a
  pure phase by this argument.
- The combination $(1-m_P)\lambda + m_P$ is exactly the massive-overlap numerator that keeps the operator
  on a rescaled GW circle; that is why this particular ratio, and not a naive $\det(D+m)$, is the one that
  is a phase.

## So why does the PDF only claim it "in the weak field limit"?

Because $|\lambda-1|=1$ is exact only if $V=\gamma\,\mathrm{sgn}(H_W)$ is **exactly** unitary. In the code
the sign function is a Zolotarev rational approximation with a finite number of poles (`npole=21`). The
approximation is excellent where $H_W$ is well-conditioned (eigenvalues bounded away from 0) and degrades
near zero modes of $H_W$. So:

- The eigenvalues drift slightly **off** the GW circle by an amount set by the Zolotarev residual.
- $|R|-1$ is then $O(\text{Zolotarev error})$, smallest in weak/smooth fields (no small $H_W$ modes) and
  largest in rough/strong fields (near-zero $H_W$ modes, where $\mathrm{sgn}$ is hardest to approximate).

This is exactly the free-field number from `run_reweighting_R_claude.log`:

$$
R_\text{free} = (0.7234,\,-0.6904), \qquad |R_\text{free}| = 1 + 5.0\times 10^{-8},
$$

i.e. a phase to $\sim 5\times10^{-8}$, the residual being the npole=21 Zolotarev error, not a genuine
modulus. The PDF's hedge "at least in the weak field limit" is the honest statement that on a general
background the only thing protecting $|R|=1$ is how well the sign function (hence the GW circle) is
realized numerically.

## Predicted result of the Gaussian-field check

For a general field drawn from a Gaussian, $|R|$ should stay equal to 1 up to the Zolotarev residual,
**provided** $H_W$ has no near-zero modes. As the Gaussian width grows the fields get rougher, $H_W$
develops small eigenvalues, the GW circle is realized less accurately, and $|R|-1$ should grow and
correlate with the maximal off-circle deviation $\max_i\big||\lambda_i-1|-1\big|$. Diagnostic to report
per sample: $|R|-1$ alongside $\max_i\big||\lambda_i-1|-1\big|$ and $\min_i|\lambda_i|$ (the
would-be zero mode).
