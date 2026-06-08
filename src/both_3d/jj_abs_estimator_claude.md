# Estimating $|G(\Delta t)|$ with error bars from the stochastic estimator

How to get the magnitude of the (real) current correlator -- for log plots / power-law (CFT) fits --
without the noise-floor bias, using the existing stochastic outputs.  Companion to
`jj_reality_tsymmetry_claude.md`.

## 0. What we know about the signal (the key fact)

Per sample $i$ (a config for an interacting ensemble; a hit for the free field) the program gives a
**complex** number $C_i(\Delta t)$.  We established:

- the true correlator $G(\Delta t)$ is **real**, and the gauge-averaged imaginary part is **strictly zero**
  (exact: lattice $T$-symmetry (3.58) $+$ trace cyclicity).  So

$$
\frac1N\sum_i \mathrm{Re}\,C_i \xrightarrow{N\to\infty} G(\Delta t),\qquad
\frac1N\sum_i \mathrm{Im}\,C_i \xrightarrow{N\to\infty} 0\ \text{(no signal at all -- pure noise).}
$$

Two consequences we will use: (i) the **signal lives entirely in $\mathrm{Re}\,C$**; (ii) $\mathrm{Im}\,C$ is
a **signal-free noise channel** we get for free.

## 1. Why $|\cdot|$ is tricky (the bias)

$|\cdot|$ is non-linear, so it does NOT commute with averaging.  For a noisy estimate $x=G+\varepsilon$ with
$\langle\varepsilon\rangle=0$, $\mathrm{Var}(\varepsilon)=\sigma^2$,
$$
\langle x^2\rangle = G^2+\sigma^2 \quad\Rightarrow\quad \langle|x|\rangle > |G|\ \text{ when } |G|\lesssim\sigma .
$$
So a magnitude built from noisy inputs is biased **upward** -- this is exactly the "noise floor" that fills
the tail of $|G|$.  Two WRONG recipes that bake in this bias:

- $\frac1N\sum_i |C_i|$  (take $|\cdot|$ per sample, then average) -- worst: each sample carries full noise.
- $\big|\frac1N\sum_i \mathrm{Re}\,C_i\big|$ evaluated where the mean is consistent with $0$ -- biased by the
  mean's own variance $\sigma^2$.

## 2. The estimator (do this)

Let $r_i=\mathrm{Re}\,C_i(\Delta t)$.

**Step 1 -- average first (linear, unbiased).**  $\bar r=\frac1N\sum_i r_i$.  This estimates $G$ (signed).

**Step 2 -- jackknife.**  Block/jackknife over samples to get the central $\bar r$ and its error $\sigma$
($\sigma^2=$ jackknife variance of $\bar r$).  Keep the SIGN -- do all resampling on the signed $\bar r$.

**Step 3 -- magnitude with bias subtraction.**  Because $\langle\bar r^2\rangle=G^2+\sigma^2$,
$$
\boxed{\ \widehat{|G|}(\Delta t)=\sqrt{\max\!\big(0,\ \bar r^{\,2}-\sigma^2\big)}\ }
$$
- **signal region** ($\bar r^2\gg\sigma^2$): $\widehat{|G|}\approx|\bar r|$, and the subtraction is negligible;
- **noise floor** ($\bar r^2\lesssim\sigma^2$): the argument clips to $0$ -> quote an **upper limit**
  $|G|\lesssim\sigma$ (or $2\sigma$), NOT a value.

**Step 4 -- error bar on $\widehat{|G|}$.**
- Signal region: propagate, $\delta\widehat{|G|}\approx \sigma\,|\bar r|/\widehat{|G|}\ (\to\sigma$ for
  $\bar r^2\gg\sigma^2$).  Equivalently jackknife $|\bar r_{(j)}|$ over resamples $j$ (gives the same band
  where there is signal).
- Floor region: the meaningful output is the upper limit $\sim\sigma$, not a central value $\pm$ error.

(If you prefer a single fully-jackknifed object: on each resample $j$ form
$m_{(j)}=\sqrt{\max(0,\bar r_{(j)}^2-\sigma^2)}$ with the GLOBAL $\sigma^2$, then take the jackknife mean/error
of $\{m_{(j)}\}$.  In the signal region this agrees with the propagation; at the floor it returns $\approx0$
with $\sim\sigma$ spread, i.e. the upper limit.)

## 3. The imaginary part as a free noise gauge (it has zero signal)

Since $\langle\mathrm{Im}\,C\rangle=0$ exactly, $\mathrm{Im}\,C$ is a clean, signal-free probe:

1. **Consistency / bug check.**  $\overline{\mathrm{Im}\,C}$ must be $0$ within its jackknife error at every
   $\Delta t$.  A nonzero value beyond noise signals a bug (or, for an $m_P$ ensemble, the genuine parity
   signal -- but for massless / $m_F$ it must vanish).
2. **Noise calibration / where the floor is.**  If the stochastic source noise is isotropic in the complex
   plane ($Z_2$ / Gaussian sources -> $\mathrm{Var}(\mathrm{Im}\,C)\approx\mathrm{Var}(\mathrm{Re}\,C)$),
   the Im channel gives an INDEPENDENT estimate of the per-sample noise, hence of $\sigma$ and of the
   $\Delta t$ where $|\bar r|$ is "just noise."  Concretely: $|\bar r|\lesssim \mathrm{std}(\overline{\mathrm{Im}\,C})$
   flags the floor.
3. Do NOT fold Im into the signal -- it carries no $G$; use it only as the gauge above.

## 4. Pseudo-code (numpy / notebook)

```
# C : array (nsamp, Nt) complex stochastic estimator (signed); binsize for jackknife
re = C.real                       # signal channel
im = C.imag                       # signal-FREE noise channel

jkR = jk_of(re)                   # signed jackknife of Re  -> mean gbar(Dt), err sig(Dt)
gbar, sig = jkR.mean(), jkR.err()

absG = np.sqrt(np.maximum(0.0, gbar**2 - sig**2))     # bias-subtracted magnitude
dabsG = np.where(absG > 0, sig*np.abs(gbar)/np.maximum(absG,1e-300), sig)  # error (signal); ~sig at floor
is_floor = gbar**2 <= sig**2                          # mark these as upper limits ~sig

# noise gauge (must be ~0 for massless/m_F):
jkI = jk_of(im); imbar, imerr = jkI.mean(), jkI.err()
# floor cross-check: |gbar| <~ std over samples of Im mean  => noise-dominated
```

Plot $\widehat{|G|}$ with $\delta\widehat{|G|}$ on a log axis; draw the floor points as downward arrows
(upper limits).  Overlay the CFT power law there.

## 5. Golden rule

Average the **signed / linear** estimator first; apply $|\cdot|$, $\log$, or $m_\text{eff}$ to the resampled
**means**, never per sample.  The notebook's signed plotter already obeys this for the band; this note adds
(a) the $-\sigma^2$ debias for the central magnitude and (b) the Im-based floor/consistency gauge.
