# IT388_Group_Project: Parallel Protein Folding Simulation - Srida Kalidindi, Srinath Shankar, James Strickert

**Project Overview**:
=====================
A parallel Monte Carlo search for low-energy protein folds in the 2D HP model, written in C with OpenMP, MPI, and a hybrid MPI + OpenMP version, and benchmarked up to 128 cores on SDSC Expanse.

*The Problem*
- A protein is a string of **H** (hydrophobic) and **P** (polar) residues placed on a 2D square grid. The chain stays connected and no two residues share a square. 
- **Energy = -(number of H-H pairs on neighboring squares that are not neighbors in the chain). Lower is better.**
- The number of possible folds grows exponentially with length, so we search for good folds instead of listing them all. 

*Approach*
- One **search** is a Monte Carlo run: 
    1. Start from a straight chain. 
    2. Make a random *pull move* 
    3. Keep the new fold if its energy is lower
    4. Repeat for N **steps** and remember the best fold (lowest energy) seen.

The steps inside one search depend on each other, but **different searches are independent**. That is where the parallelism is: we run 128 searches per protein and spread them over threads and/or processes. Each search's random seed depends only on its protein and search number, never on the thread or rank, so **every version must find exactly the same best energies as the sequential version**. The job scripts check this automatically.

*How the Work is Split*
| Version    | Program        | How Work is Split                                                     | 
|------------|----------------|-----------------------------------------------------------------------|
| Sequential | bin/sequential | One core runs every search of every protein. Baseline for speedup.|
| OpenMP     | bin/openmp     | Proteins one at a time; the searches of each protein are split across threads.|
| MPI        | bin/mpi        | Proteins are dealt round robin. Each rank runs all searches for its proteins.|
| Hybrid     | bin/hybrid     | MPI deals out the proteins, inside each rank OpenMP splits the searches. Cores = ranks x threads.|

**Contribution Breakdown:**
===========================

*Srinath*: hybrid.c , mpi.c , folding.h
*James*: folding.c , openmp.c , sequential.c , folding.h
*Srida*: proteins_known.txt , proteins.txt , collect_results.sh , 4 jobscripts
*Group*: final presentation , README.md

**Initial File Structure & Purpose**:
====================================
```text
IT388_Group_Project/
├── README.md
├── NOTES.md
├── src/
│   ├── folding.h
│   ├── folding.c
│   ├── sequential.c
│   ├── openmp.c
│   ├── mpi.c
│   └── hybrid.c
├── data/
│   ├── proteins.txt
│   └── proteins_known.txt
└── scripts/
    ├── pfold_test_known_seq.sb
    ├── pfold_test_known_par.sb
    ├── pfold_seq.sb
    ├── pfold_par.sb
    └── collect_results.sh
```
*Purpose*
| File                               | Purpose        | 
|------------------------------------|----------------|
| src/folding.h                      | Shared data types (Protein, FoldResult) and function list. Fixed-size arrays only, so a protein can be sent with one MPI broadcast. |
| src/folding.c                      | The folding core used by every version: one Monte Carlo search, the pull move, the energy calculation, the fold validity check, and the data file reader. |
| src/sequential.c                   | Baseline. All searches one after another; prints the results table and timings. |
| src/openmp.c                       | Splits the searches of each protein across threads (#pragma omp parallel for). |
| src/mpi.c                          | Rank 0 reads the file and broadcasts the proteins; ranks fold their proteins; results are combined on rank 0 with MPI_Reduce. |
| src/hybrid.c                       | mpi.c plus an OpenMP loop over the searches inside each rank. The final version. |
| data/proteins.txt                  | All 64 proteins, sorted by length. |
| data/proteins_known.txt            | 8 benchmark proteins (S20 to S64) with published energies, for the fold-quality check. |
| scripts/pfold_test_known_seq.sb    | Quick sequential test on the known proteins. |
| scripts/pfold_test_known_par.sb    | Quick parallel test on the known proteins. |
| scripts/pfold_seq.sb               | The full sequential experiment on the 64 proteins. |
| scripts/pfold_par.sb               | The full parallel experiment on the 64 proteins. |
| scripts/collect_results.sh         | Combines and performs calculations on all run files into two CSV files. |

*Data File structure*
| File | Proteins | Used for |
|---|---|---|
| `proteins.txt` | 64: 11 benchmarks + 53 random | Correctness tests and all speed/scaling runs |
| `proteins_known.txt` | 8 benchmarks (S20 to S64) | Fold-quality check against published energies |

*Data Format*: one protein per line, three space-separated fields.
```
name   sequence               best_known_energy
S20    HPHPPHHPHPPHPHHPPHPH   -9
R20_1  PHHPHHHHHPPHHPHPHPPH   0
```
- **name:** a label, up to 31 characters, no spaces.
- **sequence:** the chain, written as `H` (hydrophobic) and `P` (polar) letters, up to 127 residues.
- **best_known_energy:** the lowest energy published for that sequence. Lower (more negative) is better. The programs only print it next to the energy they found (`Best_E` vs `Known_E`). It is never used in the calculation.

*Two kinds of proteins, told apart by the first letter of the name:*
- **`S...` (11 benchmarks):** standard test sequences from ref[1] with published best energies. The number is the chain length (`S20` has 20 residues); `S100A` and `S100B` are two different 100-residue sequences.
- **`R...` (53 random):** sequences we generated for a bigger workload (random H/P, about 50% H, lengths 20-100, fixed seed so they are reproducible). In `R57_28`, 57 is the length and 28 is just a counter. They have no published energy, so their third column is **0, meaning "unknown"**, not a real energy. Only the `S` proteins can be used to judge fold quality.

*Rules:* the file is sorted by length (so the MPI round-robin split spreads the work evenly), at most 64 proteins are read, and **no comment or blank header lines are allowed**: the parser treats every line as a protein.

**Running the Program**:
========================

The experiments are run on **Expanse**, from the **project folder.**
run the following commands in this order:

1. *Setup*
module purge; module load slurm cpu/0.17.3b gcc/10.2.0 openmpi/4.1.3
mkdir -p bin results/raw results/known_test results/slurm
gcc   -g -Wall -O2 -o bin/sequential src/sequential.c src/folding.c -lm
gcc   -g -Wall -O2 -fopenmp -o bin/openmp src/openmp.c src/folding.c -lm
mpicc -g -Wall -O2 -o bin/mpi src/mpi.c src/folding.c -lm
mpicc -g -Wall -O2 -fopenmp -o bin/hybrid src/hybrid.c src/folding.c -lm
ls bin        # should list: hybrid  mpi  openmp  sequential

2. *Test the Known Proteins - run from the project folder*
sbatch scripts/pfold_test_known_seq.sb
sbatch scripts/pfold_test_known_par.sb
cat results/slurm/pfold_test_known_par.*.txt     # want: PASS for openmp, mpi and hybrid. This tells us the parallel version found exactly the same best energy as the sequential version for every protein, and every fold passed the validity check.

3. *Run the Full Experiment on the 64 Proteins*
sbatch scripts/pfold_seq.sb
sbatch scripts/pfold_par.sb

4. *After Both Jobs Finish, Build the Summary Files*
bash scripts/collect_results.sh

*The Job Scripts:*
| Script | Runs | Resources/Time Limit | Output Folder |
|--------|------|----------------------|---------------|
|pfold_test_known_seq.sb | Sequential on proteins_known.txt, 32 searches x 20,000 steps | 1 core, 10 min | results/known_test/ |
|pfold_test_known_par.sb | Sequential reference, then OpenMP (4 threads), MPI (4 ranks), hybrid (2 ranks x 4 threads), then PASS/FAIL | 8 cores, 15 min | results/known_test/ |
|pfold_seq.sb | Sequential on proteins.txt: 128 searches, 10,000 and 50,000 steps, 3 repeats (6 runs) | 1 core, 1 hour | results/raw/ |
|pfold_par.sb | OpenMP (1-128 threads), MPI (1-128 ranks), hybrid (11 ranks x threads combinations up to 128 cores), same steps and repeats (162 runs) | 1 full node (128 cores), 2.5 hours | results/raw/ |

**Output Files Structure & Purpose**:
=====================================
```text
results/
├── slurm/           job logs (what each job printed)
├── known_test/      quick tests on the known proteins
├── raw/             one file per run of the full experiment, plus checks
├── summary.csv      one row per run
└── summary_avg.csv  one row per configuration
```
*One Run File*: 
>results/raw/<version>_r<ranks>_t<threads>_n<steps>_s<searches>_rep<repeat>.out

    - name says what was run: openmp_r1_t16_n50000_s128_rep2.out = OpenMP, 1 rank, 16 threads, 50,000 steps, 128 searches, repeat 2. Sequential files start with seq_
    - A table with one row per protein: Name, Len, Best_E (best energy found), Known_E (published best; 0 = unknown), Valid (fold passed the check). MPI and hybrid add Rank (which rank folded it). Time(s) is that protein's time.

*Other Files:*
| File | What it is |
| --- | --- |
| results/raw/energy_check.txt | Every parallel run's best energies compared with the sequential run of the same size, e.g. `162 runs match, 0 mismatch, 0 not compared yet`. Differences are listed as `MISMATCH`. "Not compared yet" only means the sequential job had not finished; whichever job finishes last redoes the check. |
| results/raw/ENV_<jobid>.txt | Date, computer name and CPU (for the experimental setup). |
| results/raw/failures.log | Exists only if a run failed verification. |
| results/known_test/known_seq.out | Sequential run on the known proteins: `Best_E` next to `Known_E`, the quality check. |
| results/known_test/known_openmp_t4.out, known_mpi_r4.out, known_hybrid_r2_t4.out | The parallel test runs. `known_par_ref_seq.out` is the sequential reference they are compared with. |
| results/known_test/known_*.txt | Trimmed copies (name, best energy, valid) used only for the PASS/FAIL comparison. |
| results/slurm/<jobname>.<jobid>.<node>.txt | Everything the job printed; check the newest with `ls -t results/slurm \| head`. |

*Summary Files:*
Made by `scripts/collect_results.sh`.

| File | What it is |
| --- | --- |
| summary.csv | One row per run (168 in the default sweep). Can be used to see how much repeats differ. |
| summary_avg.csv | One row per configuration, averaged over the repeats. Used to perform graphical analysis using listed calculations. |

*Purpose of summary_avg.csv Columns:*
| Column | Meaning |
| --- | --- |
| version, ranks, threads, cores, steps, searches, reps | The configuration. `cores` = ranks x threads. |
| mean_total_s, min_total_s, max_total_s | Wall-clock time over the repeats. |
| speedup | Mean sequential time (same steps and searches) / `mean_total_s`. |
| efficiency | `speedup` / `cores` (1.0 = perfect scaling). |
| mean_moves_per_sec | Monte Carlo steps per second = proteins x searches x steps / time. |
| mean_imbalance, mean_comm_s | MPI and hybrid only: load imbalance and communication/waiting time. |
| pct_moves, pct_energy | Share of the available core time (time x cores) spent making moves and computing energy; the rest is everything else, including waiting and idle cores. |
| all_verified | 1 if every repeat passed verification. |

*summary.csv Columns:*
`summary.csv` has the per-run columns:
`version, ranks, threads, cores, proteins, searches, steps, rep, total_s, moves_s, energy_s, comm_s, imbalance, verified, moves_per_sec, file`.

**References**:
===============
[1] C. Thachuk, A. Shmygelska, and H. H. Hoos, “A replica exchange Monte Carlo algorithm for
    protein folding in the HP model,” BMC Bioinformatics, vol. 8, art. 342, Sep. 2007.
    doi:10.1186/1471-2105-8-342

[2] Bowman, Gregory. “Protein Folding .” Folding@Home, foldingathome.org/science/folding.
    Accessed 15 Oct. 2026.

[3] Message Passing Interface Forum, “MPI: A Message-Passing Interface Standard, Version
    4.1,” Nov. 2023. https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report.pdf

[4] OpenMP Architecture Review Board, “OpenMP Application Programming Interface,
    Version 6.0,” Nov. 2024. https://www.openmp.org/wp-content/uploads/OpenMP-
    API-Specification-6-0.pdf

**AI Usage Acknowledgements**:
=============================
In line with the project's AI use scale, generative AI tools were used at various levels throughout the project.
Specific AI usages and their corresponding levels are documented in the header comments of the respective files.