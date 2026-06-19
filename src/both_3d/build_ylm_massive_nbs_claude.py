#!/usr/bin/env python3
# Generate one ylm-tower analysis notebook PER MASSIVE (Real-m / m_F) mass, from the massless
# template jj_ylm_prod_analysis_claude.ipynb (which is HAND-EDITED -- never overwritten here).
#
# Per the user: "remove the commented cells, but otherwise keep the exact format."
#   - drop every fully-commented / empty CODE cell (disabled cells)
#   - keep all other cells verbatim (source, id, structure)
#   - repoint ONLY the mass-dependent strings: the esn() dir glob (sea-mass ensemble dir + valence
#     suffix) and the "(gsq=8, massless)" title labels -> "$m_F=<m>$"; also fix the stale tb8 doc
#     text to tb2 (the code already globs corr_ylm_disc_tb2_nhits1).
#   - clear outputs (massless plots) so each notebook is a fresh template for its own mass.
#
# ONE-TIME generator (like build_validate_nbs_claude.py): once a per-mass notebook is hand-edited,
# do NOT re-run this -- it would overwrite those edits.
import json
import copy

SRC = "jj_ylm_prod_analysis_claude.ipynb"

# (dir-string mRe/vmRe value , title label)
MASSES = [
    ("0.200000", "0.2"),
    ("0.100000", "0.1"),
    ("0.050000", "0.05"),
    ("0.010000", "0.01"),
]

def is_comment_cell(c):
    if c["cell_type"] != "code":
        return False
    nonempty = [l for l in c["source"] if l.strip()]
    if not nonempty:
        return True
    return all(l.lstrip().startswith("#") for l in nonempty)

nb0 = json.load(open(SRC))

# cells kept once (commented/empty code cells stripped)
kept = [c for c in nb0["cells"] if not is_comment_cell(c)]
print(f"template: {len(nb0['cells'])} cells -> {len(kept)} kept ({len(nb0['cells'])-len(kept)} commented/empty dropped)")

for m_dir, m_lab in MASSES:
    nb = copy.deepcopy(nb0)
    cells = [copy.deepcopy(c) for c in kept]
    for c in cells:
        src = "".join(c["source"])
        # 1. esn(): sea-mass ensemble dir + valence suffix (longest substring first)
        src = src.replace(
            "nu01.000000nt128L1_vmRe0.000000vmIm0.000000",
            f"nu01.000000mRe{m_dir}mIm0.000000nt128L1_vmRe{m_dir}vmIm0.000000",
        )
        # 2. any remaining bare valence suffix (e.g. the prose path in the title cell)
        src = src.replace("vmRe0.000000vmIm0.000000", f"vmRe{m_dir}vmIm0.000000")
        # 3. title / heading label
        src = src.replace("massless", f"$m_F={m_lab}$")
        # 4. stale disc-dir doc text (code already uses tb2)
        src = src.replace("tb8", "tb2")
        c["source"] = src.splitlines(keepends=True)
        # clear outputs -> fresh template
        if c["cell_type"] == "code":
            c["outputs"] = []
            c["execution_count"] = None
    nb["cells"] = cells
    out = f"jj_ylm_prod_analysis_mRe{m_lab}_claude.ipynb"
    json.dump(nb, open(out, "w"), indent=1)
    print(f"wrote {out}  ({len(cells)} cells, mRe={m_dir})")
