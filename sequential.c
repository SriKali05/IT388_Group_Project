/*
 * ---------------------------------------------------------------
 *  Name:      James Strickert
 *  Course:    IT 388 - Parallel Processing
 *  Project:   Parallel Protein Folding Simulation
 *  File:      sequential.c
 *
 *  Purpose:   Baseline version. Reads the proteins, runs #searches
 *             Monte Carlo searches of #steps each on every protein,
 *             one after another, and prints the best energy and the
 *             time. This is the reference for speedup and for the
 *             quality check.
 *
 *  Compile:   gcc -g -Wall -O2 -o sequential src/sequential.c src/folding.c -lm
 *  Run:       ./sequential <#searches> <#steps> [protein file]
 *             (protein file defaults to data/proteins.txt)
 * ---------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "folding.h"

#define MAX_PROTEINS 64

double getTime(void);

int main(int argc, char *argv[])
{
    Protein proteins[MAX_PROTEINS];
    FoldResult result, best;
    int num_proteins, num_searches, num_steps;
    int i, s, all_valid = 1;
    double start, protein_start, elapsed;
    double time_moves = 0.0, time_energy = 0.0;
    const char *data_file = "data/proteins.txt";

    //Read and check the input
    if (argc != 3 && argc != 4)
    {
        fprintf(stderr, "\n---USAGE: ./sequential <#searches> <#steps> [protein file]---\n\n");
        exit(1);
    }
    num_searches = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    if (argc == 4)
        data_file = argv[3];
    if (num_searches <= 0 || num_steps <= 0)
    {
        fprintf(stderr, "\n--- #searches and #steps must be positive integers ---\n\n");
        exit(1);
    }

    num_proteins = parse_proteins(data_file, proteins, MAX_PROTEINS);
    if (num_proteins == 0)
    {
        fprintf(stderr, "No proteins loaded.\n");
        exit(1);
    }

    printf("Sequential: %d proteins, %d searches x %d steps each\n\n",
           num_proteins, num_searches, num_steps);
    printf("%-6s %4s %7s %8s %6s %10s\n", "Name", "Len", "Best_E", "Known_E", "Valid", "Time(s)");
    printf("----------------------------------------------\n");

    start = getTime();
    for (i = 0; i < num_proteins; i++)
    {
        protein_start = getTime();

        //run every search; the seed depends only on the protein and the
        //search number, so every version gets the same results
        for (s = 0; s < num_searches; s++)
        {
            result = fold_protein(&proteins[i], num_steps, (unsigned int)(i * num_searches + s + 1));
            time_moves += result.time_moves;
            time_energy += result.time_energy;
            if (s == 0 || result.best_energy < best.best_energy)
                best = result;
        }

        //quality check: the best fold must be valid and its energy must match
        int valid = check_fold(&proteins[i], best.coords) &&
                    compute_energy(&proteins[i], best.coords) == best.best_energy;
        if (!valid)
            all_valid = 0;

        printf("%-6s %4d %7d %8d %6s %10.4f\n", proteins[i].name, proteins[i].len,
               best.best_energy, proteins[i].best_known_energy,
               valid ? "yes" : "NO", getTime() - protein_start);
    }
    elapsed = getTime() - start;

    printf("\nTotal time: %f sec\n", elapsed);
    printf("Time in moves: %f sec, in energy: %f sec, other: %f sec\n",
           time_moves, time_energy, elapsed - time_moves - time_energy);
    if (all_valid)
        printf("Verification: PASSED (every best fold is valid)\n");
    else
        printf("Verification: FAILED (some best fold is not valid)\n");
    return 0;
}

//wall clock time in seconds
double getTime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
