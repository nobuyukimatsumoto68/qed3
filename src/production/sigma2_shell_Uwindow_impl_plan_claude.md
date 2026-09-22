# Shell kernel from the distillation eigenvectors $U$ (index windows), NOT from $M$

## Goal (NM, 2026-09-19)

The shell single-meson kernel must be built from $U$ = eigenvectors of the symmetrized hermitian
operator $\tilde D^\dagger\tilde D+\tilde D\tilde D^\dagger$ (the `BASIS_SYM` distillation basis), NOT
from $M=\tau(t,t)-\tfrac12$. Kernel for a window $w$ of eigenvector indices:
$$
K_w(t)=\sum_{i\in w}u_i(t)\,u_i(t)^\dagger .
$$
Since the perambulator is already in this basis, the mode-space vertex is the 0/1 diagonal block
$E_w$, and
$$
C_{ab}(t,s)=-\,\mathrm{Re}\sum_{i\in a}\sum_{j\in b}\tau(t,s)_{ij}\,\tau(s,t)_{ji}.
$$
The SLICING pattern (which index windows) is a separate input (`WINDOWS` env); it is set by the
eigenvalue clusters of the symmetrized operator itself, independent of $M$.

Check requested: at L1 (sym perambulators exist, `distill_Nv24_sym`), does "the same problem"
(shell diagonals sinking to $m_{PS}$ / contaminated $(2,2)$) occur with this kernel?

## Files

1. `sigma2_shell_Uwindow_claude.py` (NEW) -- shell-only driver: index-window kernel, diagonal
   effmasses + shell GEVP (fixed-$t_0$, rebased, no Hankel). `NVDIR` selects the basis
   (`distill_Nv24_sym` vs `distill_Nv24` for comparison). Asserts `/evals` ascending.

## Chunks

- Chunk 1 (Files: `sigma2_shell_Uwindow_claude.py`): driver + L1 run on the sym basis.
- Chunk 2: same driver on the standard basis for A/B; later extend to the combined basis
  (cross with $\sigma^2$ via $K_w=U E_w U^\dagger$) once the slicing is fixed.

## Open question

- Slicing for the symmetrized operator: its FREE degeneracy pattern is not on disk (`data_free`
  holds only the one-sided $\tilde D^\dagger\tilde D$ basis, free degeneracies 8,4,12). The
  interacting ensemble-averaged sym spectrum is smooth (pairs, no sharp cluster edges). Need either
  a free `BASIS_SYM` run (GPU, NM runs) or NM's slicing choice.

## Refs

- Distillation: Peardon et al. arXiv:0905.2160. Physics context: Chester, Pufu arXiv:1603.05582.
