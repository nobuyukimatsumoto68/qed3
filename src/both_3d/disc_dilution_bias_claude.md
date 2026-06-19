# Disconnected-loop two-point: time-dilution self-contamination (derivation + fix)

## 0. TL;DR

Time dilution makes the **single-current loop** $J(t)$ unbiased (your iid argument is exactly right).
It also makes the **two-point** $\langle J(t_0)J(t_0+\Delta t)\rangle$ unbiased for **almost all** $\Delta t$.
It fails **only** at $\Delta t$ equal to a multiple of the dilution period $P=N_t/t_\text{block}$ (here $P=16$),
and there the bias is a full **connected** propagator term that is huge compared to the true tiny disconnected
signal -- the observed spike at $\Delta t=16$. It is **not** a copied-random-numbers bug; the noise is iid in
$t$ (verified, Sec. 5). The contamination comes from using **one source / one solve to supply both legs**
at the two timeslices that happen to live in the same dilution class. Fix = make the two legs use
**independent** noise (Sec. 6).

## 1. What we want

Per gauge config, the disconnected loop at timeslice $t$ is the trace
$$
J(t) = \mathrm{Tr}_t\!\big[\Gamma\, D^{-1}\big] = \sum_{i\in t}\big(\Gamma D^{-1}\big)_{ii},
$$
where $i=(t,x,\alpha)$ runs over spatial site $x$ and spin $\alpha$ on timeslice $t$, and $\Gamma$ is the
projected current insertion ($\sigma_a$ with the $A_n Y_{\ell m}$ weight). The disconnected current-current
correlator is the (gauge-ensemble) correlation of two such loops,
$$
C(\Delta t) = \frac{1}{N_t}\sum_{t_0} J(t_0)\,J(t_0+\Delta t),
$$
then averaged over configurations. $J(t)$ is a deterministic number once the config is fixed; the only
issue is the **stochastic estimation** of $J(t)$ at every $t$.

## 2. Stochastic estimator and time dilution (what the code does)

Complex $Z_2$ noise $\eta$ with
$$
E[\eta_i]=0,\qquad E[\eta_i^{*}\eta_j]=\delta_{ij},\qquad E[\eta_i\eta_j]=0 .
$$
Solve once $\phi = D^{-1}\eta$. The loop estimator at timeslice $t$ reads the noise back **only on $t$**:
$$
\hat J(t) = \sum_{i\in t}\eta_i^{*}\,(\Gamma\phi)_i
          = \sum_{i\in t}\sum_{k}\eta_i^{*}\,\Gamma_{ij}\,(D^{-1})_{jk}\,\eta_k .
$$
(This is exactly `accumulate_loop_raw`: `J[t] += sum_x conj(eta(t,x,spin)) * (Gamma phi)(t,x,spin)`.)

**Time dilution.** Partition the $N_t$ timeslices into $P=N_t/t_\text{block}$ classes; class $c$ is the set
$$
\mathcal{T}_c=\{\,t : t\equiv c \ (\mathrm{mod}\ P)\,\}=\{c,\ c+P,\ c+2P,\dots\},\qquad c=0,\dots,P-1 .
$$
For class $c$ a **fresh independent** noise vector $\eta^{(c)}$ is drawn, supported on $\mathcal{T}_c$, and a
**separate solve** $\phi^{(c)}=D^{-1}\eta^{(c)}$ is done. The estimate $\hat J(t)$ for $t\in\mathcal{T}_c$ uses
$\eta^{(c)},\phi^{(c)}$. Crucially:
$$
c\neq c' \;\Rightarrow\; \eta^{(c)} \text{ and } \eta^{(c')} \text{ independent}.
$$
In the code `t_block=8` $\Rightarrow P=16$: each class has 8 timeslices spaced by 16, one solve per class.

## 3. The single loop is unbiased (your argument)

$$
E[\hat J(t)] = \sum_{i\in t}\sum_k \Gamma_{ij}(D^{-1})_{jk}\,E[\eta_i^{*}\eta_k]
            = \sum_{i\in t}\Gamma_{ij}(D^{-1})_{ji}
            = \mathrm{Tr}_t[\Gamma D^{-1}] = J(t).
$$
Only the $k=i$ term survives because $E[\eta_i^{*}\eta_k]=\delta_{ik}$. All the off-diagonal noise
(other timeslices in the same class) has zero mean. So summing several timeslices into one solve to save
inversions is fine **for the trace**. This is the standard dilution gain, and it is correct.

## 4. The two-point: clean almost everywhere, biased at $\Delta t=\,$multiple of $P$

We estimate
$$
\hat C(\Delta t)=\frac{1}{N_t}\sum_{t_0}\hat J(t_0)\,\hat J(t_0+\Delta t).
$$
Let $t_1=t_0+\Delta t$, with $t_0\in\mathcal{T}_{c_0}$, $t_1\in\mathcal{T}_{c_1}$,
$c_0=t_0\bmod P$, $c_1=t_1\bmod P$.

### Case A: $\Delta t \not\equiv 0 \pmod P$  (the generic case, e.g. $\Delta t=1,\dots,15$)

Then $c_0\neq c_1$, so $\hat J(t_0)$ and $\hat J(t_1)$ are built from **independent** noise vectors
$\eta^{(c_0)},\eta^{(c_1)}$. The expectation factorizes:
$$
E[\hat J(t_0)\hat J(t_1)] = E[\hat J(t_0)]\,E[\hat J(t_1)] = J(t_0)\,J(t_1).
$$
**Unbiased.** This is precisely your "cross-timeslice expectation vanishes by iid": the spurious term
would need $E[\eta^{(c_0)}\eta^{(c_1)}]$, which is zero. Good.

### Case B: $\Delta t \equiv 0 \pmod P$  (e.g. $\Delta t=16,32,\dots$)

Now $c_0=c_1\equiv c$: **both legs use the same $\eta^{(c)}$ and the same $\phi^{(c)}$.** Write the product
with all indices in class $c$:
$$
\hat J(t_0)\hat J(t_1)=
\Big(\sum_{i\in t_0}\eta_i^{*}\Gamma_{ij}D^{-1}_{jk}\eta_k\Big)
\Big(\sum_{l\in t_1}\eta_l^{*}\Gamma_{lm}D^{-1}_{mn}\eta_n\Big),
$$
with $k,n$ summed over **all** class-$c$ DOF (where $\eta$ is nonzero). Taking $E$ over the four noise
factors $\eta_i^{*},\eta_k,\eta_l^{*},\eta_n$, the two leading contractions are:

- **(i) Disconnected pairing** $\;\langle\eta_i^{*}\eta_k\rangle\langle\eta_l^{*}\eta_n\rangle=\delta_{ik}\delta_{ln}$:
$$
\sum_{i\in t_0}\Gamma_{ij}D^{-1}_{ji}\;\sum_{l\in t_1}\Gamma_{lm}D^{-1}_{ml}
= J(t_0)\,J(t_1)\qquad\text{(the signal we want).}
$$

- **(ii) Cross / "connected" pairing** $\;\langle\eta_i^{*}\eta_n\rangle\langle\eta_l^{*}\eta_k\rangle=\delta_{in}\delta_{lk}$:
this forces $n=i\in t_0$ and $k=l\in t_1$. Because **both $t_0$ and $t_1$ belong to class $c$**, those
$\eta$ components are nonzero, and the term survives:
$$
\boxed{\;\sum_{i\in t_0}\sum_{l\in t_1}\Gamma_{ij}D^{-1}_{jl}\,\Gamma_{lm}D^{-1}_{mi}
= \mathrm{Tr}\!\big[\Gamma\,D^{-1}_{t_0,t_1}\,\Gamma\,D^{-1}_{t_1,t_0}\big]\;\neq\;0.\;}
$$

So
$$
E[\hat C(\Delta t)] = C(\Delta t) \;+\; \underbrace{\frac{1}{N_t}\sum_{t_0}\mathrm{Tr}\big[\Gamma D^{-1}_{t_0,t_0+\Delta t}\Gamma D^{-1}_{t_0+\Delta t,t_0}\big]}_{\text{bias, only at }\Delta t\equiv0\ (\mathrm{mod}\ P)} .
$$

In Case A this cross pairing needs $\delta_{in}$ linking the two **different** sources, i.e.
$\langle\eta^{(c_0)}\eta^{(c_1)}\rangle=0$ -- so it is absent. The bias exists **only** when the two legs
share a source, i.e. at $\Delta t$ a multiple of $P$.

## 5. It is not a copied-noise bug

`rng.h`: `sites` is an `[Nt][n_sites]` array of independent `SingleRng`s (line 100), each reseeded
independently (lines 119-120); `CZ2_site(t,ix)` draws from `sites[t][ix]` (line 115). Hence
$\eta(t)$ and $\eta(t+P)$ are **independent draws**, not copies. The coupling between timeslices is purely
through the single solve $\phi=D^{-1}\eta$ (the $D^{-1}_{t_0,t_1}$ in term (ii)), not through identical
random numbers. So nothing in `time_spin_dilution` needs a "regenerate" fix -- the noise is already iid in $t$.

## 6. Why it is large, and the fix

The bias (ii) is a genuine **connected** propagator round trip $t_0\!\to\!t_1\!\to\!t_0$. It is the same
order as the connected correlator and therefore **enormous** compared to the true disconnected signal
(which is $\langle J\rangle\langle J\rangle$-type, tiny). That is why $\Delta t=16$ sticks out by $\sim 30\times$.

The cross term (ii) vanishes **iff the two legs use independent noise**. Options, with cost relative to the
current disc ($P$-class dilution, $n_h=1$):

1. **Cross-hit estimator (recommended, cheapest correct).** Run $n_h\ge 2$ independent hits (each with the
   same cheap time dilution). Build the two-point from **different** hits only:
   $$
   \hat C(\Delta t)=\frac{1}{N_t}\frac{1}{n_h(n_h-1)}\sum_{t_0}\sum_{h\neq h'}\hat J^{(h)}(t_0)\,\hat J^{(h')}(t_0+\Delta t).
   $$
   Since $\langle\eta^{(h)}\eta^{(h')}\rangle=0$ for $h\neq h'$, term (ii) is gone at **all** $\Delta t$.
   This is exactly the $[\,S(t_0)S(t_1)-Q(t_0,t_1)\,]/(n_h(n_h{-}1))$ "drop the $h{=}h'$ diagonal" formula
   already documented (commented) in `jj_disc_postproc_claude.cc:69-71`. Cost: $\times n_h$ (e.g. $\times2$).

2. **Full time dilution** ($t_\text{block}=1\Rightarrow P=N_t$). Each timeslice is its own class/solve, so for
   every $\Delta t\neq0$ the two legs are independent -> no bias anywhere. Cost: $N_t$ solves/config
   ($\times P = \times16$ here). Exact but expensive.

3. **Mask** $\Delta t\in\{P,2P,\dots\}$. Zero cost, but leaves holes; not a real fix.

Reference for the unbiased disconnected estimator (drop the diagonal hit): standard stochastic-source
technique, e.g. Bali, Collins, Schafer, Comput. Phys. Commun. 181 (2010) 1570 [arXiv:0910.3970];
McNeile-Michael, Phys. Rev. D 73 (2006) 074506.

## 7. Recommendation

Use **option 1**: bump the disc to $n_h\ge 2$ and form the two-point from distinct hits (drop the diagonal).
It keeps the cheap dilution, removes the bias at every $\Delta t$, and matches the postproc's intended
(currently disabled) estimator. $n_h=2$ already kills the bias; more hits just lower the variance.
The connected (conn) estimator is unaffected -- it is a genuine source-to-sink propagator and has no such
self-contraction.
