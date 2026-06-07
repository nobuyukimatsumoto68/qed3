# `jj_corr_claude.cu` -- clean specification (equation <-> code)

The **combined** current-current correlator program: connected **and** disconnected in one pass, sharing
$\phi'=D_m^{-1}\eta$ and the $K(\cdot,t)\phi'$ sink applies.  Output: `data_<ESNID>/corr_nt0<N>_nhits<H>/corr.<k>.h5`.

This file is `jj_conn_correlators_claude.cu` (connected only) **plus the disconnected fold** described in
Sec. A below.  The connected estimators -- operators, weights, shared forward leg, vector tp/sp/ylm tower,
axial $C_{A+-}$ -- are **identical** to `jj_conn_correlators_spec_claude.md` Sections 1-5; only the output
directory differs (`corr_` instead of `conn_`).  Read that spec first; this one specifies the delta.

Reference: `qed3int_v2-10.pdf`; design history `conserved_current_correlators_impl_plan_v3_claude.md` Sec. 3.8.

---

## A. Disconnected single-current traces (the fold)

The disconnected diagram is the projected single-current expectation
$$
J(t)=\sum_a w_a\,T(a,t),\qquad T(a,t)=\eta^\dagger K(a,t)\,\phi',\quad \phi'=D_m^{-1}\eta,
$$
with $E[T(a,t)]=\mathrm{tr}\big(K(a,t)\,D_m^{-1}\big)$.  **Raw** (no $1/4\pi$ fold; `write_vec` `:356`),
matching the standalone `jj_disc_claude.cu`.  Origin-independent: accumulated over the whole hit, written once.

**Reuse (Sec. 3.8 lever).**  $T(a,t)$ needs exactly the sink quantity $K(a,t)\phi'$ that the connected $(++)$
sink passes already compute.  So the disc traces RIDE those loops at **zero extra $K$ applies**:

| projection | accumulator | rides | code |
|---|---|---|---|
| tp  | $J^{tp}(t)=\sum_n w^{tp}_n\,\eta^\dagger\,\texttt{kphi}$ | tp+ylm sink (`kphi`$=K(n,t)\phi'$) | `:445` |
| ylm | $J_\ell(t)=\eta^\dagger\,\Phi_\ell(t)$, $\Phi_\ell=\sum_n W_\ell(n)\texttt{kphi}$ | same loop (after $\Phi_\ell$ accumulation) | `:449` |
| sp  | $J^{sp}(t)=\sum_a w^{sp}_{il}\,\eta^\dagger\,\texttt{kphi}$ | sp sink (`kphi`$=K(\ell_a,t)\phi'$) | `:526` |

Accumulators declared per hit (`:411`); written at `:570-572` as `h{h}/disc/{tp,sp,ylm/l{l}}/J`.

**Parity tilde trace.**  For $m$ imaginary the dagger-leg single-trace is
$$
\tilde T(a,t)=\big(K(a,t)\,\tilde\phi\big)^\dagger\eta,\qquad \tilde\phi=\tilde D_{m_P}^{-1}\eta .
$$
This CANNOT ride the connected parity sink (which applies $K^\dagger\,\texttt{phimm}$ with
$\texttt{phimm}=\tilde D_{m_P}^{-\dagger}\eta$) -- different operator direction and vector.  So it gets its own
forward solve + sink pass (`:574-601`): `tilphi` solve (`:579`), then $K(\cdot,t)\tilde\phi$ over sites/links,
written as `h{h}/disc/{tp,sp,ylm/l{l}}/Jtil`.  (Massless / $m_F$: $\tilde T=T^*$, reconstructed downstream,
not dumped.)

---

## B. Output (h5; `corr.<k>.h5` in `data_<ESNID>/corr_nt0<N>_nhits<H>/`)

Connected keys exactly as `jj_conn_correlators_spec_claude.md` Sec. 6 (under the `corr_` dir).  Added per hit:
```
h{h}/disc/tp/J/{real,imag}          h{h}/disc/sp/J/{real,imag}
h{h}/disc/ylm/l{l}/J/{real,imag}                                          l=0,1,2
h{h}/disc/tp/Jtil/{real,imag}       h{h}/disc/sp/Jtil/{real,imag}         (parity only)
h{h}/disc/ylm/l{l}/Jtil/{real,imag}                                       (parity only)
```
Disc keys are RAW (no $1/4\pi$) and length $N_t$, origin-independent (no `t0_{b}` level).
Resume gate unchanged: the `"complete"` sentinel (written LAST) -- read-only check before recompute.

---

## C. Three-program map (run script `run_jj_analysis_claude.sh`)

| flag | binary | output dir | content |
|---|---|---|---|
| `--all` (default) | `jj_corr_claude.o` | `corr_nt0<N>_nhits<H>/` | connected + folded disc (this spec) |
| `--conn` | `jj_conn_correlators_claude.o` | `conn_nt0<N>_nhits<H>/` | connected only |
| `--disc` | `jj_disc_claude.o` | `disc_nhits<H>/` | disconnected only (cheap; no connected solves) |

`--all` is the efficient production path (one $\phi'$, shared sink).  `--conn`/`--disc` are the standalone
single-diagram programs, kept for separate / cheap-many-config runs.

---

## D. Solve count (per hit) vs. the connected-only program

Connected source solves are unchanged from `jj_conn_correlators_spec_claude.md` Sec. 7.  The fold adds:
**0** solves in the non-parity case (disc rides the existing sink applies); **1** solve in the parity case
(`tilphi`) plus its $(n_\text{sites}+n_\text{links})N_t$ tilde-sink $K$ applies.
