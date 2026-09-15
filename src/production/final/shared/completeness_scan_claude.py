# completeness_scan_claude.py
# COMMON-STYLE checkpoint-count (ncfg) table generator, channel-agnostic.
# Reference implementation (qed3-d4 / Conn A); extended to cover glue (qed3-a7, 2026-09-14).
#
# Counts per-config h5 files for the standard massless grid, applying the kmin=20 thermalization
# cut.  Glob only -- no h5 reads.  Files may live in a per-config SUBDIR under the ensemble dir
# (fermion channels) or DIRECTLY in the ensemble dir (glue).
#
# Usage examples:
#   axial (conn only, fermion valence dir + corr.<k>.h0.h5):
#     python3 completeness_scan_claude.py --conn-subdir corr_ylm_conn_t00_nhits1_s1
#   vector (conn + disc):
#     python3 completeness_scan_claude.py --conn-subdir corr_ylm_conn_t00_nhits1_s1 \
#         --disc-subdir corr_ylm_disc_tb2_nhits1
#   glue (bare gauge dir, files directly in the ensemble dir, glue_f2_v2_shapes.<k>.h5):
#     python3 completeness_scan_claude.py --bare --conn-subdir "" \
#         --file-prefix glue_f2_v2_shapes --file-suffix .h5 --label glue
#
# Common conventions (propagate verbatim):
#   - kmin = 20 thermalization cut (discard config index k < 20).  SEPARATE from the ncfg>=100
#     minimum post-cut count used for the LOW flag.
#   - columns: L | Nf | gsq | at | ncfg_conn | ncfg_disc | status
#     (a conn-only channel prints "n/a" in ncfg_disc; --disc-subdir omitted).
#   - status keyed on the RELEVANT count(s): OK (>=100) / LOW(<100) / EMPTY / NO_CONNDIR / NO_ENSDIR.
#   - ordering: L ascending, then gsq ascending, then Nf ascending.
#
# Channel-shape knobs (defaults reproduce the original fermion behavior byte-for-byte):
#   --bare          use the BARE gauge dir data_Nf..._hb<hb>/ (no _vmRe.._vmIm.. valence suffix).
#                   Fermion channels omit this (valence _vm dir is the default).
#   --conn-subdir   per-config subdir under the ensemble dir; pass "" when files sit DIRECTLY in the
#                   ensemble dir (glue).  Default "" -- fermion callers pass their subdir explicitly.
#   --file-prefix / --file-suffix  per-config file pattern <prefix>.<k><suffix>
#                   (default prefix "corr", suffix ".h0.h5" -> corr.<k>.h0.h5; glue: "glue_f2_v2_shapes"/".h5").
import argparse
import glob
import os
import re

Nt = 128
KMIN = 20
LOW_STATS_NCFG = 100
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
NFS = [2, 4, 6]


def esn(nf, L, g, at, bare):
    base = ("data_Nf%d_gsq%.6fat%.6fnu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            % (nf, g, at, L, HB[L]))
    if bare:
        return base + "/"
    return base + "_vmRe0.000000vmIm0.000000/"


def count_kmin(searchdir, prefix, suffix, kmin):
    fs = glob.glob(searchdir + "/" + prefix + ".*" + suffix)
    n = 0
    for f in fs:
        m = re.search(re.escape(prefix) + r"\.(\d+)\.", f)
        if m is None:
            continue
        if int(m.group(1)) >= kmin:
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="common-style ncfg completeness table")
    ap.add_argument("--conn-subdir", default="",
                    help="per-config conn subdir under the ensemble dir; \"\" = files directly in the ensemble dir")
    ap.add_argument("--disc-subdir", default="", help="per-config disc subdir (omit for conn-only channels)")
    ap.add_argument("--bare", action="store_true",
                    help="use the bare gauge dir (no _vm valence suffix); default off = fermion valence dir")
    ap.add_argument("--file-prefix", default="corr", help="per-config file prefix (glue: glue_f2_v2_shapes)")
    ap.add_argument("--file-suffix", default=".h0.h5", help="per-config file suffix (glue: .h5)")
    ap.add_argument("--at", type=float, default=0.2)
    ap.add_argument("--kmin", type=int, default=KMIN)
    ap.add_argument("--label", default="channel", help="channel label for the header")
    args = ap.parse_args()

    have_disc = len(args.disc_subdir) > 0
    print("# %s ensemble completeness (at=%.1f, massless) -- kmin=%d thermalization cut"
          % (args.label, args.at, args.kmin))
    print("# ncfg = per-config h5 with index k >= kmin. status: OK(>=%d)/LOW(<%d)/EMPTY/NO_CONNDIR/NO_ENSDIR."
          % (LOW_STATS_NCFG, LOW_STATS_NCFG))
    print()
    print("| L | Nf | gsq | at | ncfg_conn | ncfg_disc | status |")
    print("|---|----|-----|----|-----------|-----------|--------|")
    nrows = 0
    nlow = 0
    for L in [1, 2, 3, 4]:
        for g in GS[L]:
            for nf in NFS:
                ensdir = esn(nf, L, g, args.at, args.bare)
                cdir = ensdir + args.conn_subdir
                ncfg_conn = count_kmin(cdir, args.file_prefix, args.file_suffix, args.kmin) if os.path.isdir(cdir) else 0
                if have_disc:
                    ddir = ensdir + args.disc_subdir
                    ncfg_disc = count_kmin(ddir, args.file_prefix, args.file_suffix, args.kmin) if os.path.isdir(ddir) else 0
                    disc_s = str(ncfg_disc)
                    ref = min(ncfg_conn, ncfg_disc)   # matched analysis is conn-disc limited
                else:
                    disc_s = "n/a"
                    ref = ncfg_conn
                if not os.path.isdir(ensdir):
                    status = "NO_ENSDIR"
                elif not os.path.isdir(cdir):
                    status = "NO_CONNDIR"
                elif ref == 0:
                    status = "EMPTY"
                elif ref < LOW_STATS_NCFG:
                    status = "LOW(<%d)" % LOW_STATS_NCFG
                    nlow += 1
                else:
                    status = "OK"
                print("| %d | %d | %.1f | %.1f | %d | %s | %s |"
                      % (L, nf, g, args.at, ncfg_conn, disc_s, status))
                nrows += 1
    print()
    print("# roll-up: %d ensembles listed, %d LOW(<%d)." % (nrows, nlow, LOW_STATS_NCFG))


if __name__ == "__main__":
    main()
