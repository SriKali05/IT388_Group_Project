/*
 * ---------------------------------------------------------------
 *  Name:      James Strickert
 *  Course:    IT 388 - Parallel Processing
 *  Project:   Parallel Protein Folding Simulation
 *  File:      openmp.c
 *
 *  Purpose:   OpenMP version. Same work as sequential.c, but the
 *             independent searches on one protein are split across
 *             threads. Each search uses the same seed as in the
 *             sequential version, so the best energies must match it
 *             exactly.
 *
 *  Compile:   gcc -g -Wall -O2 -fopenmp -o openmp src/openmp.c src/folding.c -lm
 *  Run:       ./openmp <#threads> <#searches> <#steps> [protein file]
 *             (protein file defaults to data/proteins.txt)
 * 
 *  AI Usage ACK: In line with the project's AI use scale, Claude, a generative AI tool, was used at a level 4 (co-writer) to turn our 
 *  parallel algorithm design into pseudocode. 
 * ---------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "folding.h"

#define MAX_PROTEINS 64

int main(int argc, char *argv[])
{
    Protein proteins[MAX_PROTEINS];
    FoldResult best;
    int thread_count, num_proteins, num_searches, num_steps;
    int i, s, best_search, all_valid = 1;
    double start, protein_start, elapsed;
    double time_moves = 0.0, time_energy = 0.0;
    const char *data_file = "data/proteins.txt";

    //Read and check the input
    if (argc != 4 && argc != 5)
    {
        fprintf(stderr, "\n---USAGE: ./openmp <#threads> <#searches> <#steps> [protein file]---\n\n");
        exit(1);
    }
    thread_count = atoi(argv[1]);
    num_searches = atoi(argv[2]);
    num_steps = atoi(argv[3]);
    if (argc == 5)
        data_file = argv[4];
    if (thread_count <= 0 || num_searches <= 0 || num_steps <= 0)
    {
        fprintf(stderr, "\n--- #threads, #searches and #steps must be positive integers ---\n\n");
        exit(1);
    }

    num_proteins = parse_proteins(data_file, proteins, MAX_PROTEINS);
    if (num_proteins == 0)
    {
        fprintf(stderr, "No proteins loaded.\n");
        exit(1);
    }

    printf("OpenMP: %d proteins, %d searches x %d steps each, threads: %d\n\n",
           num_proteins, num_searches, num_steps, thread_count);
    printf("%-6s %4s %7s %8s %6s %10s\n", "Name", "Len", "Best_E", "Known_E", "Valid", "Time(s)");
    printf("----------------------------------------------\n");

    start = omp_get_wtime();
    for (i = 0; i < num_proteins; i++)
    {
        protein_start = omp_get_wtime();
        best_search = -1;

        //searches are independent, so each thread runs its share of them
#pragma omp parallel for num_threads(thread_count) reduction(+:time_moves, time_energy)
        for (s = 0; s < num_searches; s++)
        {
            //same seed as sequential.c: depends only on the protein and search number
            FoldResult result = fold_protein(&proteins[i], num_steps, (unsigned int)(i * num_searches + s + 1));
            time_moves += result.time_moves;
            time_energy += result.time_energy;

            //keep the lowest energy; on a tie keep the lower search number,
            //which is the one the sequential version keeps
#pragma omp critical
            {
                if (best_search == -1 || result.best_energy < best.best_energy ||
                    (result.best_energy == best.best_energy && s < best_search))
                {
                    best = result;
                    best_search = s;
                }
            }
        }

        //quality check: the best fold must be valid and its energy must match
        int valid = check_fold(&proteins[i], best.coords) &&
                    compute_energy(&proteins[i], best.coords) == best.best_energy;
        if (!valid)
            all_valid = 0;

        printf("%-6s %4d %7d %8d %6s %10.4f\n", proteins[i].name, proteins[i].len,
               best.best_energy, proteins[i].best_known_energy,
               valid ? "yes" : "NO", omp_get_wtime() - protein_start);
    }
    elapsed = omp_get_wtime() - start;

    printf("\nTotal time: %f sec\n", elapsed);
    printf("Time in moves: %f sec, in energy: %f sec (summed over all %d threads)\n",
           time_moves, time_energy, thread_count);
    if (all_valid)
        printf("Verification: PASSED (every best fold is valid)\n");
    else
        printf("Verification: FAILED (some best fold is not valid)\n");
    return 0;
}
