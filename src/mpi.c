/*
 * ---------------------------------------------------------------
 *  Name:      Srinath Shankar
 *  Course:    IT 388 - Parallel Processing
 *  Project:   Parallel Protein Folding Simulation
 *  File:      mpi.c
 *
 *  Purpose:   MPI version. Rank 0 reads the proteins and sends them
 *             to every rank. Proteins are dealt out round robin
 *             (protein i goes to rank i % nprocs) so the long proteins
 *             at the end of the file don't all end up on the last rank.
 *             Each rank runs all the searches for its own proteins one
 *             after another, then the results are combined on rank 0.
 *             Seeds are the same as in sequential.c, so the best
 *             energies must match the sequential version exactly.
 *
 *  Compile:   mpicc -g -Wall -O2 -o mpi src/mpi.c src/folding.c -lm
 *  Run:       mpirun -n <#procs> ./mpi <#searches> <#steps> [protein file]
 *             (protein file defaults to data/proteins.txt)
 * 
 *  AI Usage ACK: In line with the project's AI use scale, Claude, a generative AI tool, was used at a level 4 (co-writer) to turn our 
 *  parallel algorithm design into pseudocode. 
 * ---------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "folding.h"

#define MAX_PROTEINS 64

int main(int argc, char *argv[])
{
    Protein proteins[MAX_PROTEINS];
    FoldResult result, best;
    int rank, nprocs;
    int num_proteins = 0, num_searches, num_steps;
    int i, s, all_valid = 1;
    const char *data_file = "data/proteins.txt";

    //results for the proteins this rank owns (the rest keep the start values)
    int my_energy[MAX_PROTEINS], my_valid[MAX_PROTEINS];
    double my_ptime[MAX_PROTEINS];
    double my_breakdown[2] = {0.0, 0.0};    //time in moves, time in energy

    //combined results, only used on rank 0
    int all_energy[MAX_PROTEINS], all_valid_list[MAX_PROTEINS];
    double all_ptime[MAX_PROTEINS];
    double total_breakdown[2];
    double *rank_times = NULL;

    double start, elapsed, t, protein_start;
    double comp_time, comm_time = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //Read and check the input (every rank has argv, so no need to send it)
    if (argc != 3 && argc != 4)
    {
        if (rank == 0)
            fprintf(stderr, "\n---USAGE: mpirun -n <#procs> ./mpi <#searches> <#steps> [protein file]---\n\n");
        MPI_Finalize();
        exit(1);
    }
    num_searches = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    if (argc == 4)
        data_file = argv[3];
    if (num_searches <= 0 || num_steps <= 0)
    {
        if (rank == 0)
            fprintf(stderr, "\n--- #searches and #steps must be positive integers ---\n\n");
        MPI_Finalize();
        exit(1);
    }

    //only rank 0 reads the file
    if (rank == 0)
    {
        num_proteins = parse_proteins(data_file, proteins, MAX_PROTEINS);
        if (num_proteins == 0)
        {
            fprintf(stderr, "No proteins loaded.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    //start the clock once everyone is ready
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //send the proteins to every rank
    t = MPI_Wtime();
    MPI_Bcast(&num_proteins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(proteins, num_proteins * (int)sizeof(Protein), MPI_BYTE, 0, MPI_COMM_WORLD);
    comm_time += MPI_Wtime() - t;

    //start values: energy 1 is higher than any real energy so MPI_MIN skips it,
    //valid 1 so MPI_MIN only drops to 0 if the owner found a bad fold
    for (i = 0; i < num_proteins; i++)
    {
        my_energy[i] = 1;
        my_valid[i] = 1;
        my_ptime[i] = 0.0;
    }

    //fold my proteins: rank, rank + nprocs, rank + 2*nprocs, ...
    t = MPI_Wtime();
    for (i = rank; i < num_proteins; i += nprocs)
    {
        protein_start = MPI_Wtime();

        //same seed as sequential.c: depends only on the protein and search number
        for (s = 0; s < num_searches; s++)
        {
            result = fold_protein(&proteins[i], num_steps, (unsigned int)(i * num_searches + s + 1));
            my_breakdown[0] += result.time_moves;
            my_breakdown[1] += result.time_energy;
            if (s == 0 || result.best_energy < best.best_energy)
                best = result;
        }

        //quality check is done here since only this rank has the coords
        my_energy[i] = best.best_energy;
        my_valid[i] = check_fold(&proteins[i], best.coords) &&
                      compute_energy(&proteins[i], best.coords) == best.best_energy;
        my_ptime[i] = MPI_Wtime() - protein_start;
    }
    comp_time = MPI_Wtime() - t;

    //combine everything on rank 0
    //(this time also includes rank 0 waiting for the slowest rank to finish)
    t = MPI_Wtime();
    MPI_Reduce(my_energy, all_energy, num_proteins, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(my_valid, all_valid_list, num_proteins, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(my_ptime, all_ptime, num_proteins, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(my_breakdown, total_breakdown, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    comm_time += MPI_Wtime() - t;
    elapsed = MPI_Wtime() - start;

    //compute time of each rank, to see how even the work split is
    if (rank == 0)
        rank_times = malloc(nprocs * sizeof(double));
    MPI_Gather(&comp_time, 1, MPI_DOUBLE, rank_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("MPI: %d proteins, %d searches x %d steps each, ranks: %d\n\n",
               num_proteins, num_searches, num_steps, nprocs);
        if (nprocs > num_proteins)
            printf("Note: more ranks than proteins, %d ranks have no work\n\n", nprocs - num_proteins);

        printf("%-6s %4s %7s %8s %6s %5s %10s\n", "Name", "Len", "Best_E", "Known_E", "Valid", "Rank", "Time(s)");
        printf("----------------------------------------------------\n");
        for (i = 0; i < num_proteins; i++)
        {
            if (!all_valid_list[i])
                all_valid = 0;
            printf("%-6s %4d %7d %8d %6s %5d %10.4f\n", proteins[i].name, proteins[i].len,
                   all_energy[i], proteins[i].best_known_energy,
                   all_valid_list[i] ? "yes" : "NO", i % nprocs, all_ptime[i]);
        }

        printf("\nTotal time: %f sec\n", elapsed);
        printf("Time in moves: %f sec, in energy: %f sec (summed over all %d ranks)\n",
               total_breakdown[0], total_breakdown[1], nprocs);
        printf("Communication time on rank 0: %f sec (includes waiting for other ranks)\n", comm_time);

        //load balance: slowest rank vs the average rank
        double max_time = 0.0, sum_time = 0.0;
        printf("Compute time per rank:");
        for (i = 0; i < nprocs; i++)
        {
            printf(" %.4f", rank_times[i]);
            sum_time += rank_times[i];
            if (rank_times[i] > max_time)
                max_time = rank_times[i];
        }
        printf("\nLoad imbalance (slowest / average): %.2f\n", max_time / (sum_time / nprocs));

        if (all_valid)
            printf("Verification: PASSED (every best fold is valid)\n");
        else
            printf("Verification: FAILED (some best fold is not valid)\n");
        free(rank_times);
    }

    MPI_Finalize();
    return 0;
}
