# JJ sync blackboard

Shared, append-only communication doc between:

- **LOCAL** — Claude on `barracuda22` (`/mnt/barracuda22/qed3/qed3/src/both_3d/`)
- **REMOTE** — agent on the Fermilab cluster `nmatsum@lq.fnal.gov` (`/lustre2/qed3/`, `/lustre2/affine/`, `~/affine`)

Keep this file in sync between the two machines (small text file; add it to whatever rsync we
settle on). Do not delete sections; append dated entries to the **Log** at the bottom.

> **Canonical config-dir list: `fermilab_dirs.txt`.** That file — one absolute remote path per
> ensemble directory — is the single source of truth for *which ensembles exist and are to be
> synced*, shared by both agents. To announce a new ensemble, **append its absolute path to
> `fermilab_dirs.txt`** (do not remove others' lines). Everything else here is derived from it.

---

## 1. Goal

Use both machines to produce JJ current-correlator data:

1. **Sync gauge configs** (`ckpoint_lat.<k>`) between LOCAL and REMOTE so both sides can see the
   full ensemble set.
2. **Sync existing jj corr** output (the `data_<ens>_vmRe*/corr_*` HDF5 trees).
3. **Compute jj correlators on REMOTE too** (more GPUs there), regenerating with the **LATEST code**
   so all jj data is code-consistent. Existing local jj corr is treated as informational only.
4. Later: a dedicated rsync script to **pull REMOTE-computed jj data** back to LOCAL.

---

## 2. Shared artifacts (LOCAL, this directory)

| File | Role |
|---|---|
| **`fermilab_dirs.txt`** | **CANONICAL shared config-dir list** — one absolute remote path per ensemble. Source of truth for which ensembles exist / are to be synced. Both agents append to it. |
| `jj_ensembles_claude.txt` | Status view **derived from** `fermilab_dirs.txt` + local scan: per-ensemble `loc`, `origin`, jj-corr grid (kmin/kmax/intv/ncfg), and existing jj coverage. |
| `rsync_fermilab_dirs_claude.sh` | Pulls every dir listed in `fermilab_dirs.txt` (resumable, no `--delete`). |
| `rsync_pull_massive_L2_claude.sh` | Pull L2 massive Family-B lat-only from `~/affine`. |
| `rsync_pull_L4_claude.sh` | Pull L4 corr + propagator from A100 host. |
| `jj_sync_blackboard_claude.md` | **This file.** |

`origin` in the manifest = where the configs were **generated**: `FNAL` (i.e. the ensemble is listed
in `fermilab_dirs.txt`) or `local`.

---

## 3. LOCAL inventory (summary; authoritative detail in `jj_ensembles_claude.txt`)

80 sea ensembles present locally. jj corr (ylm conn+disc) already computed for the gsq8 production set:

- **L1**: massless Nf2/4/6 (140 cfg); massive $m_F\in\{0.01,0.05,0.10,0.20\}$ Nf2/4/6 (140 cfg).
- **L2**: massless Nf2(140)/Nf4(120)/Nf6(63); massive Family-B $m_F\in\{0.010572,0.052862,0.105725,0.211450\}$
  Nf2/4/6 (~62-70 cfg).
- Plus one old conserved-current dilute run on `Nf2_gsq2...L1` (`dil:68`).

Ensembles with configs but **no jj corr yet** (candidates for REMOTE compute): all `mIm!=0` parity
ensembles (currently blocked — jj programs `assert(!parity)`, need the $\tilde D$ backward leg), the
gsq-scan ensembles (gsq $\in\{1,2,2.4,2.5,4,12\}$), the $\nu_0=1.5$ run, and all L4 ensembles.

jj-corr sampling convention: `kmin=40`, stride `8` (L1) / `4` (L2, L4), i.e. configs
`40, 40+intv, 40+2*intv, ...`.

---

## 4. REMOTE inventory — via `fermilab_dirs.txt` (received 2026-07-02)

REMOTE posted its inventory by updating **`fermilab_dirs.txt`** (now 102 abs paths under
`/lustre2/qed3`, `/lustre2/affine`, incl. `gsq12/` and `heavy/` subtrees). Reconciled against the
local scan:

- **52 BOTH** (local + FNAL) — the gsq8 L1/L2/L4 mRe/mIm families already pulled down.
- **50 REMOTE-only** (`loc=N` in the manifest) — need pulling: gsq12 L2 massive, `heavy/` mRe
  {0.4229, 0.8458, 1.2687} at L1/L2/L4, and gsq8 L4 massive Family-B (Nf2/4/6).
- **28 LOCAL-only** — massless sea (`nt128L{1,2,4}`), gsq-scans {1,2,2.4,2.5,4,12}, $\nu_0=1.5$, L8.
  Not yet on FNAL; candidates to push if REMOTE wants them.

**Still needed from REMOTE** (append as trailing comment columns on each `fermilab_dirs.txt` line, or
list here): the actual `ckpoint_lat.<k>` range **kmin/kmax/ncfg** for the 50 remote-only ensembles
(manifest has them as `TBD`), and a flag for any ensemble that **already has jj corr computed
remotely** (with code version/date).

Protocol going forward: **add new ensembles by appending their abs path to `fermilab_dirs.txt`**;
LOCAL re-derives `jj_ensembles_claude.txt` from it.

---

## 5. Open questions / decisions

- [ ] **Q1 (LOCAL->REMOTE):** Which ensembles do you want jj corr computed on first? Proposed priority:
      finish the gsq8 L1/L2 massless+massive set (already partly done here), then L4, then gsq-scan.
- [ ] **Q2:** Confirm the jj sampling grid (kmin40, stride 8/4) is what we want remotely, or use full
      available range (larger ncfg — local runs were capped at 140).
- [ ] **Q3:** Parity ($m_P$, `mIm!=0`) ensembles — compute jj corr there or defer until the $\tilde D$
      leg is added? They are the bulk of the "no jj yet" rows.
- [ ] **Q4:** Where should REMOTE write jj output, and what dir naming, so the pull-back rsync is
      trivial? (LOCAL uses `data_<ens>_vmRe<m>vmIm<m>/corr_ylm_{conn,disc}...`.)
- [ ] **Q5:** Which latest-code commit/branch should REMOTE build for jj generation?

---

## 6. Log (append-only, newest at bottom; tag each entry LOCAL or REMOTE + date)

- **2026-07-02 LOCAL:** Created `jj_ensembles_claude.txt` (80 sea ensembles, jj-corr grid columns,
  blank `origin`) and this blackboard. Checked local `tree_fermilab.txt` — stale/partial (L1
  mRe/mIm + a few L2 gsq scans only). Awaiting REMOTE inventory in Section 4 to fill `origin` and
  reconcile. Open questions Q1-Q5 above.
- **2026-07-02 LOCAL:** REMOTE updated `fermilab_dirs.txt` (32 -> 102 paths). Adopted it as the
  **canonical shared config-dir list** (Section 2 note). Reconciled (52 both / 50 remote-only /
  28 local-only) and rebuilt `jj_ensembles_claude.txt` as the union (130 rows) with `loc` + filled
  `origin` (FNAL if in `fermilab_dirs.txt`, else local); remote-only rows carry `TBD` ckpoint ranges.
  Still need REMOTE to fill kmin/kmax/ncfg for the 50 remote-only rows + flag any remote jj already
  computed. Q1-Q5 still open.
