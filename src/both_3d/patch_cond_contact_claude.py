#!/usr/bin/env python3
# SURGICAL one-off (NOT the builder).  Subtracts the contact term in the two condensate plot cells of the
# three per-mass notebooks (condensate_contact_massive_claude.md): contact = 2 (forward 1/2 + backward 1/2,
# x tr[A1]/Vst=2), with the backward leg rescaled to 1/(1+m) -> contact = 1+1/(1+m) for sigma_PS in the m_P
# run and sigma_FS in the m_F run.  Done via TARGETED line replacements (not whole-cell) so any hand-edits
# (e.g. the user's `print(e/Vst)`) and outputs are preserved.  Guarded + idempotent (skips if 'csub' present).
import json

def patch(src, sig, eq, rescale):
    if 'csub' in src:
        return src, False                     # already contact-subtracted
    tex = '\\' + sig.replace('_', '_{') + '}'  # \sigma_{PS}
    R = [
        ("# CONDENSATE %s: stoch (Re,Im) vs dense exact.  /Vst." % sig,
         "# CONDENSATE %s, CONTACT-SUBTRACTED (condensate_contact_massive_claude.md): sigma/Vst + |contact|." % sig),
        ("e,er,ei = jk_scalar(raw['%s']); dd = dense['%s']/Vst" % (sig, sig),
         "e,er,ei = jk_scalar(raw['%s'])\n"
         "m_cur = {'massless':0.0,'mF':0.1+0j,'mP':0.1j}[MASS]; csub = (1.0+1.0/(1.0+m_cur)) if MASS=='%s' else 2.0  # |contact|\n"
         "dd = dense['%s']/Vst + csub" % (sig, rescale, sig)),
        ("plt.errorbar([0,1],[e.real/Vst,e.imag/Vst],yerr=[er/Vst,ei/Vst],fmt='o',capsize=4,label='stoch (jk)')",
         "es = e/Vst + csub\n"
         "plt.errorbar([0,1],[es.real,es.imag],yerr=[er/Vst,ei/Vst],fmt='o',capsize=4,label='stoch (jk)')"),
        ("plt.xticks([0,1],['Re','Im']); plt.ylabel(r'$\\langle%s\\rangle / V_{st}$')" % tex,
         "plt.axhline(0,color='gray',lw=0.5); plt.xticks([0,1],['Re','Im']); plt.ylabel(r'$\\langle%s\\rangle/V_{st}-$contact')" % tex),
        ("r'$%s$ (%s): stoch vs dense, $m='" % (tex, eq),
         "r'$%s$ (%s) contact-subtracted, $m='" % (tex, eq)),
        ("/cond_%s_'+MLAB" % sig,
         "/cond_%s_sub_'+MLAB" % sig),
    ]
    ok = all(a in src for a, _ in R)
    if not ok:
        return src, False                     # something hand-edited -> skip (do not clobber)
    for a, b in R:
        src = src.replace(a, b)
    return src, True

for fn in ['jj_validate_m0_claude.ipynb', 'jj_validate_mF_claude.ipynb', 'jj_validate_mP_claude.ipynb']:
    nb = json.load(open(fn)); n = 0; skipped = []
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        src = c['source'] if isinstance(c['source'], str) else ''.join(c['source'])
        if 'CONDENSATE sigma_PS' in src:
            new, ch = patch(src, 'sigma_PS', '1.23', 'mP')
        elif 'CONDENSATE sigma_FS' in src:
            new, ch = patch(src, 'sigma_FS', '1.55', 'mF')
        else:
            continue
        if ch:
            c['source'] = new.splitlines(keepends=True); n += 1
        else:
            skipped.append(c.get('id'))
    if n:
        json.dump(nb, open(fn, 'w'), indent=1)
    print('%s: patched %d, skipped %s' % (fn, n, skipped or 'none'))
