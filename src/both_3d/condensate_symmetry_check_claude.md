# Step 1: are $\sigma_{PS}$ parity-symmetric and $\sigma_{FS}$ flavor-symmetric?

Source: `main.tex` Sec. "Parity" (Eqs. eq:parity_Psi, eq:parity_Psibar, eq:parity_symm_mass) and
Sec. "Flavor symmetry / Lattice" (Eqs. eq:def_psi, eq:def_psibar, and the $\sigma_F$ definition).
Conventions: single pair ($N_f=2$), spinor indices contracted, spatial label summed (densities; a
bilinear $\sum_x \sigma(x)$). $\xi,\eta$ are the two 2-component blocks of $\Psi=(\xi,\eta)^T$,
$\bar\Psi=(\eta^\dagger,\xi^\dagger)$.

## A. $\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta=\bar\Psi\Psi$ is PARITY-symmetric

Parity (eq:parity_Psi / eq:parity_Psibar), with $\sigma_1^2=1$:
$$
\xi(x)\to\sigma_1\eta(x_P),\qquad \eta(x)\to\sigma_1\xi(x_P),
$$
$$
\eta^\dagger(x)\to\xi^\dagger(x_P)\sigma_1,\qquad \xi^\dagger(x)\to\eta^\dagger(x_P)\sigma_1.
$$
Apply term by term:
$$
\eta^\dagger(x)\,\xi(x)\;\to\;\big[\xi^\dagger(x_P)\sigma_1\big]\big[\sigma_1\eta(x_P)\big]
   =\xi^\dagger(x_P)\,\eta(x_P),
$$
$$
\xi^\dagger(x)\,\eta(x)\;\to\;\big[\eta^\dagger(x_P)\sigma_1\big]\big[\sigma_1\xi(x_P)\big]
   =\eta^\dagger(x_P)\,\xi(x_P).
$$
So the two terms SWAP, and the density maps to itself at the reflected point:
$$
\sigma_{PS}(x)\;\to\;\xi^\dagger(x_P)\eta(x_P)+\eta^\dagger(x_P)\xi(x_P)=\sigma_{PS}(x_P).
$$
Summing over all $x$ (relabel $x\leftrightarrow x_P$, a bijection): $\sum_x\sigma_{PS}\to\sum_x\sigma_{PS}$.
=> parity-invariant. The $\sigma_1$ factors are exactly what is needed; they cancel pairwise. This is
the [Appelquist-Pisarski-Narayanan] parity mass $\bar\Psi\Psi$ (eq:parity_symm_mass). CONFIRMED.

## B. $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$ is FLAVOR-symmetric

On the lattice the flavor multiplet is (eq:def_psi / eq:def_psibar), $V\equiv D_{ov}-1$ (unitary):
$$
\psi_1=\tfrac{1}{\sqrt2}(\xi-V^\dagger\eta),\quad \psi_2=\tfrac{1}{\sqrt2}(\xi+V^\dagger\eta),\quad
\bar\psi_1=\tfrac{1}{\sqrt2}(\eta^\dagger-\xi^\dagger),\quad \bar\psi_2=\tfrac{1}{\sqrt2}(\eta^\dagger+\xi^\dagger).
$$
The flavor-singlet diagonal mass is, by direct expansion,
$$
\sigma_{F}\equiv\bar\psi_1\psi_1+\bar\psi_2\psi_2
=\tfrac12\big[(\eta^\dagger-\xi^\dagger)(\xi-V^\dagger\eta)+(\eta^\dagger+\xi^\dagger)(\xi+V^\dagger\eta)\big]
=\eta^\dagger\xi+\xi^\dagger V^\dagger\eta .
$$
Using $V^\dagger=D_{ov}^\dagger-1=-(1-D_{ov}^\dagger)$:
$$
\boxed{\;\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta\;}
$$
which is exactly the coded `etadag_xi` $-$ `xidag_1mDdag_eta` and the paper's $\sigma_F$.

Why it is flavor-symmetric: it is the SINGLET $\sum_i\bar\psi_i\psi_i$ in the multiplet basis. The flavor
$SU(2)$ acts (eq:gamma4psi/eq:gamma5psi/eq:gamma45psi) as the vector rotations
$1\otimes\tau_3,\,-1\otimes\tau_2,\,1\otimes\tau_1$, i.e. $\delta\psi=\tau\psi$, $\delta\bar\psi=-\bar\psi\tau$.
Hence
$$
\delta\big(\bar\psi_i\psi_i\big)=(-\bar\psi\tau)\psi+\bar\psi(\tau\psi)=0
$$
for each generator -> $\sigma_{FS}$ invariant. The $(1-D_{ov}^\dagger)$ is precisely the GW dressing
($V^\dagger$) that makes the singlet structure hold on the lattice; in the continuum $V\to-1$ and
$\sigma_{FS}\to\eta^\dagger\xi-\xi^\dagger\eta=\sigma_F^{\rm cont}=\bar\Psi\gamma_{4,5}\Psi$. CONFIRMED.

## Complementarity (the reason for the names, for the NEXT step)
In the SAME variables, $\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta$ is the flavor NON-singlet
($\propto\bar\chi\tau_3\chi$, continuum), and $\sigma_{FS}$ is parity-odd. So:
- $\sigma_{PS}$: parity-EVEN, flavor-breaking.
- $\sigma_{FS}$: flavor-EVEN, parity-breaking.
They are the two complementary mass terms. (To verify the "breaks the other" half explicitly is the
natural Step 2.)
