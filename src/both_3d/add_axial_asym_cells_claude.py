#!/usr/bin/env python3
# Append an AXIAL t-asymmetry diagnostic section (1 markdown + 2 code cells: tp, sp) to each per-mass
# ylm analysis notebook.  SURGICAL: existing cells (and their outputs) are preserved; we only append.
# Idempotent: skips a notebook that already has the section.
#
# Diagnostic: A(dt) = g^A(dt) - g^A(Nt-dt) for the axial m-summed Vpp tower (l=1).  Plots |Re A|, |Im A|
# on log vs (i) their own jackknife noise floor and (ii) the symmetric signal |Re g|.  |A| pinned to the
# noise floor -> consistent with zero (t-symmetric); |A| detaching above it -> genuine t-asymmetry
# (parity-odd / Chern-Simons).  Uses the notebook's NF knob + existing _load/jackknife/conn_files/Nt/at.
import json

NB_MASS = {
    "jj_ylm_prod_analysis_mRe0.2_claude.ipynb": "0.2",
    "jj_ylm_prod_analysis_mRe0.1_claude.ipynb": "0.1",
    "jj_ylm_prod_analysis_mRe0.05_claude.ipynb": "0.05",
    "jj_ylm_prod_analysis_mRe0.01_claude.ipynb": "0.01",
}

SENTINEL = "AXIAL t-asymmetry diagnostic"

def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {}, "execution_count": None,
            "outputs": [], "source": src.splitlines(keepends=True)}

def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {}, "source": src.splitlines(keepends=True)}

MD = (
    "## " + SENTINEL + "\n"
    "$A(\\Delta t)=g^A(\\Delta t)-g^A(N_t-\\Delta t)$ for the axial $\\ell=1$ tower (m-summed `Vpp`), Re and Im.\n"
    "Plotted as $|A|$ (log) against its own jackknife **noise floor** and the symmetric signal "
    "$|\\mathrm{Re}\\,g|$.  $|A|$ pinned to the noise floor $\\Rightarrow$ consistent with zero "
    "(t-symmetric); $|A|$ detaching above it $\\Rightarrow$ genuine t-asymmetry (parity-odd / "
    "Chern-Simons).  Uses the `NF` knob; re-run as stats grow."
)

# tp = s3 ; sp = s1 + s2.  CHAN is the python expression that builds g (nconf,Nt) from cf.
def body(chan_label, g_expr, mlab):
    return (
        f"# AXIAL {chan_label} t-asymmetry: |A(dt)|, A = g(dt)-g(Nt-dt), Re and Im, vs noise floor + signal.\n"
        "fig, ax = plt.subplots()\n"
        "l = 1\n"
        "dts = np.arange(1, Nt//2)\n"
        "t = dts*at\n"
        "cf = conn_files(NF)\n"
        f"g = {g_expr}\n"
        "asym = g - g[:, (Nt - np.arange(Nt)) % Nt]\n"
        "am, are_, aim_ = jackknife(asym)\n"
        "gm, gre, gim = jackknife(g)\n"
        "ax.plot(t, np.abs(gm.real[dts]), color='0.6', lw=0.9, label=r'signal $|\\mathrm{Re}\\,g|$')\n"
        "ax.plot(t, are_[dts], color='C0', lw=0.7, ls=':', label='Re noise floor')\n"
        "ax.plot(t, aim_[dts], color='C1', lw=0.7, ls=':', label='Im noise floor')\n"
        "ax.scatter(t, np.abs(am.real[dts]), s=12, color='C0', label=r'$|\\mathrm{Re}\\,A|$')\n"
        "ax.scatter(t, np.abs(am.imag[dts]), s=12, color='C1', marker='s', label=r'$|\\mathrm{Im}\\,A|$')\n"
        "ax.set_yscale('log')\n"
        "ax.set_xlim(0, at*Nt/2)\n"
        "ax.set_xlabel(r'$t=a_t n_t$')\n"
        f"ax.set_ylabel(r'$|g^{{A,{chan_label}}}(\\Delta t)-g^{{A,{chan_label}}}(N_t-\\Delta t)|$')\n"
        f"ax.set_title(r'AXIAL {chan_label} t-asymmetry ($\\ell=1$, $m_F={mlab}$), Nf='+str(NF))\n"
        "ax.legend(fontsize=7)\n"
        "plt.tight_layout()\n"
    )

TP_EXPR = "sum(_load(cf, f'h0/ylm_axial/s3/l{l}/m{m}/Vpp') for m in range(-l, l+1))/(2*l+1)"
SP_EXPR = "sum(_load(cf, f'h0/ylm_axial/s1/l{l}/m{m}/Vpp') + _load(cf, f'h0/ylm_axial/s2/l{l}/m{m}/Vpp') for m in range(-l, l+1))/(2*l+1)"

for fn, mlab in NB_MASS.items():
    nb = json.load(open(fn))
    if any(SENTINEL in "".join(c["source"]) for c in nb["cells"]):
        print(f"{fn}: already has the section -- skipped")
        continue
    nb["cells"].append(md_cell("axasym_md", MD))
    nb["cells"].append(code_cell("axasym_tp", body("tp", TP_EXPR, mlab)))
    nb["cells"].append(code_cell("axasym_sp", body("sp", SP_EXPR, mlab)))
    json.dump(nb, open(fn, "w"), indent=1)
    print(f"{fn}: appended 3 cells (now {len(nb['cells'])})")
