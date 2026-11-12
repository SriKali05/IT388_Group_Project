# NOTES

Running log of experiments, bugs, and dead ends. Newest at the bottom.
Put your name and the date on each entry.

## 2026-11-08 (James): folding core, sequential and OpenMP versions

What was built (branch `james-folding-core`, based on Srinath's `srinath-mpi-skeleton`):

- `src/folding.c`: the shared core. 2D square lattice, straight-chain start,
  pull moves (Lesh et al. 2003, the move set ref [1] uses), Metropolis rule
  at a fixed temperature (`TEMPERATURE 0.5`). One `fold_protein()` call is one
  search. Energy is the plain pair-by-pair count, O(n^2).
- `src/sequential.c`: baseline. `./sequential <#searches> <#steps> [file]`
- `src/openmp.c`: the searches of one protein split across threads.
  `./openmp <#threads> <#searches> <#steps> [file]`

Changes to the shared `folding.h` (Srinath and Srida: please read):

- `FoldResult` has two new fields, `time_moves` and `time_energy` (seconds),
  for the time breakdown. `mpi.c` still compiles without changes.
- New function `check_fold()`: returns 1 if a fold is connected and
  self-avoiding. Used for the quality check.
- `parse_proteins()` now lives in `folding.c`. **The copy in `mpi.c` has to be
  deleted**, or the MPI build fails with a duplicate definition.

Seeds: search `s` on protein `i` uses seed `i * num_searches + s + 1`. It does
not depend on the thread or rank, so every version runs exactly the same
searches and must get exactly the same best energies. `mpi.c` currently uses
`rank * 1000 + pi * 100 + s`, which changes with the number of ranks and
repeats seeds once `s >= 100`; switching to the formula above makes MPI match
too.

Checks run (WSL Ubuntu, gcc 15.2, 12 cores; prototyping only, not Expanse):

- Builds clean with `-Wall -Wextra`.
- Stress test: 6.6 million random pull moves on chains of length 2 to 64,
  under the address and undefined-behavior sanitizers. Every new fold was
  connected and self-avoiding.
- HPPH (the square example) gives energy -1.
- `openmp` on 1, 4 and 12 threads gives exactly the same best energies as
  `sequential` (100 searches x 5000 steps).

Temperature sweep, best energy with 12 threads (known: -9 -9 -8 -14 -23 -21 -36 -42):

| T | 100 x 5000 | 20 x 50000 |
|---|---|---|
| 0.25 | -9 -9 -8 -12 -15 -17 -24 -20 | -9 -9 -8 -13 -21 -17 -31 -30 |
| 0.35 | -9 -9 -8 -13 -17 -17 -27 -22 | -9 -9 -8 -14 -21 -19 -33 -31 |
| 0.5 | -9 -9 -8 -13 -19 -17 -28 -23 | -9 -9 -8 -14 -21 -21 -34 -33 |
| 0.7 | -9 -9 -7 -12 -20 -16 -30 -24 | -9 -8 -8 -13 -21 -18 -34 -35 |

- 0.5 is best overall, so it stays.
- Longer searches help much more than more searches. Short proteins reach
  the known energy at any setting; S48 to S64 need tens of thousands of steps
  per search. Worth using more steps than `mpi.c`'s `MC_STEPS 5000` for the
  real runs.

Timing:

- Sequential, 100 x 5000: 1.15 s total. About 64% energy, 24% moves,
  12% other. Energy is the hot spot.
- The per-step timers themselves cost about 15-25% (2.5 s vs 2.1 s on
  200 x 5000). Same cost in every version, so speedup stays fair, but
  absolute times are that much higher.
- OpenMP, 100 x 5000: 1.10 s (1 thread), 0.40 s (4), 0.27 s (12). This run
  is too small for real speedup numbers; those must come from Expanse.

Not done yet / ideas:

- Faster energy (only recount the residues that moved) if energy time
  limits us.
- Replica exchange (ref [1]) as the stretch goal for the long proteins.

## 2026-11-12 (Srinath): MPI and hybrid versions

- Merged main into `srinath-mpi-skeleton`. James's files were uploaded to the
  repo root, so I moved them into `src/` and dropped the root copies of
  `mpi.c` and `proteins.txt` (same as `data/proteins.txt`). Layout now matches
  the plan: `src/` for code, `data/proteins.txt` for input.
- `src/mpi.c`: deleted my copy of `parse_proteins()` (it's in `folding.c` now).
  Seeds switched to `i * num_searches + s + 1`. #searches and #steps are now
  command line args instead of #defines.
  `mpirun -n <#procs> ./mpi <#searches> <#steps> [file]`
- `src/hybrid.c`: MPI across proteins + James's OpenMP loop across searches.
  `mpirun -n <#procs> ./hybrid <#threads> <#searches> <#steps> [file]`
- Proteins go to ranks round robin (protein i -> rank i % nprocs) instead of
  blocks. The file is sorted short to long, so blocks put S60 and S64 on the
  same rank. With 4 ranks the compute time per rank was 0.10 0.10 0.15 0.18 s.
- Each rank checks its own best folds with `check_fold()` (only that rank has
  the coords). Rank 0 gets energies/valid flags with MPI_Reduce(MIN) and the
  time breakdown with MPI_Reduce(SUM), so no coords are sent.
- Both print: total time, moves/energy time, communication time on rank 0
  (includes waiting for the slowest rank), compute time per rank and a load
  imbalance number (slowest rank / average rank) for the results slides.

Checks (Ubuntu 24.04 container, gcc 13 + OpenMPI; prototyping only, not Expanse):

- All four programs build clean with `-Wall -Wextra`.
- 40 x 5000: `mpi` on 3 ranks and `hybrid` on 1x1, 4x2 and 10x2 all give
  exactly the same best energies as `sequential`, all folds valid.

Limits to talk about in the presentation:

- Only 8 proteins, so MPI can use at most 8 ranks. Extra ranks sit idle
  (the program prints a note). More cores have to come from threads.
- The run time is set by the rank that gets S64. Even round robin can't fix
  that, the longest protein is a floor on the time.
