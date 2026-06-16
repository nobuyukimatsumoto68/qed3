#!/usr/bin/env python3
# SURGICAL one-off (NOT the builder; does not regenerate).  Replaces ONLY the buggy prep_phys block in the
# disc-load cell of the three per-mass notebooks with the unified delete-1 jackknife (error now carries the
# disc variance + conn-disc covariance).  Guarded: only edits a cell that contains the EXACT old block, so
# it cannot clobber a hand-edited cell; all other cells, outputs, and commented-out cells are left untouched.
import json

OLD = (
"def prep_phys(proj):                     # physical = -C_conn + C_disc (Eq 3.39, vector V++), Nt/2-shift'd\n"
"    conn = half_shift_avg(dil_hits(DIL, proj+'/Vpp'))\n"
"    J = disc_J(proj); Jbar = J.mean(0); Cdisc = half_shift_avg(INV4PI*two_point(Jbar, Jbar))\n"
"    m,ere,eim = jk_cplx(-conn + Cdisc[None,:])               # free: -conn(hits) + disc central (~0)\n"
"    d = -half_shift_avg(mirror(EXS_SCALE*det_raw(EXS, proj+'/Vpp')))   # -conn determ (+ disc determ = 0 free)\n"
"    return m,ere,eim,d\n"
)

NEW = (
"def prep_phys(proj):                     # physical = -C_conn + C_disc (Eq 3.39, vector V++), Nt/2-shift'd.\n"
"    # SINGLE delete-1 jackknife over hits: each pseudo-value recomputes BOTH the conn mean and the\n"
"    # (quadratic) disc two-point on the SAME leave-one-out hit set, so the error carries the disc variance\n"
"    # AND the conn-disc covariance (disc is NOT added as a fixed constant -- that undercounts the error).\n"
"    conn = half_shift_avg(dil_hits(DIL, proj+'/Vpp')); J = disc_J(proj)\n"
"    H = len(conn); Sc = conn.sum(0); SJ = J.sum(0)\n"
"    samp = []\n"
"    for i in range(H):\n"
"        Jbar_i = (SJ - J[i])/(H-1)\n"
"        samp.append(-(Sc - conn[i])/(H-1) + half_shift_avg(INV4PI*two_point(Jbar_i, Jbar_i)))\n"
"    samp = np.array(samp); est = samp.mean(0)\n"
"    ere = np.sqrt((H-1)*np.mean((samp.real-est.real)**2,0))\n"
"    eim = np.sqrt((H-1)*np.mean((samp.imag-est.imag)**2,0))\n"
"    d = -half_shift_avg(mirror(EXS_SCALE*det_raw(EXS, proj+'/Vpp')))   # -conn determ (+ disc determ = 0 free)\n"
"    return est, ere, eim, d\n"
)

for fn in ['jj_validate_m0_claude.ipynb','jj_validate_mF_claude.ipynb','jj_validate_mP_claude.ipynb']:
    nb = json.load(open(fn))
    n = 0
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        src = c['source'] if isinstance(c['source'], str) else ''.join(c['source'])
        if OLD in src:
            src = src.replace(OLD, NEW)
            c['source'] = src.splitlines(keepends=True)   # store as line list (Jupyter's format)
            n += 1
    if n:
        json.dump(nb, open(fn, 'w'), indent=1)
        print('patched %s  (%d cell)' % (fn, n))
    else:
        print('NO EXACT MATCH in %s -- left untouched (already patched or hand-edited)' % fn)
