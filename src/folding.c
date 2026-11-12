/*
 * ---------------------------------------------------------------
 *  Name:      James Strickert
 *  Course:    IT 388 - Parallel Processing
 *  Project:   Parallel Protein Folding Simulation
 *  File:      folding.c
 *
 *  Purpose:   Shared folding core used by every version of the program.
 *             A protein is a chain of H and P residues on a 2D square
 *             lattice. Residue i sits at (coords[2i], coords[2i+1]).
 *             Energy = -(number of H-H pairs that touch on the lattice
 *             but are not next to each other in the chain).
 *
 *             fold_protein() runs one Monte Carlo search:
 *               1. start from a straight chain
 *               2. make a random pull move (Lesh et al. 2003, also used
 *                  by ref [1])
 *               3. keep it if the energy is lower, or sometimes if it
 *                  is higher (Metropolis rule at a fixed temperature)
 *               4. repeat num_steps times, remembering the best fold
 *
 *             Every step depends on the one before it, so one search
 *             is sequential. The parallel versions run many searches
 *             at once, each with its own seed.
 * ---------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "folding.h"

//Metropolis temperature: higher accepts more uphill moves
#define TEMPERATURE 0.5

//the 4 lattice directions: right, up, left, down
static const int DX[4] = {1, 0, -1, 0};
static const int DY[4] = {0, 1, 0, -1};

static double getTime(void);
static int isOccupied(const int *coords, int n, int x, int y);
static int isAdjacent(int x1, int y1, int x2, int y2);
static int pullMove(const int *cur, int *trial, int n, unsigned int *seed);

//run one Monte Carlo search on protein p and return the best fold found
FoldResult fold_protein(const Protein *p, int num_steps, unsigned int seed)
{
    FoldResult result;
    int cur[MAX_SEQ_LEN * DIM];
    int trial[MAX_SEQ_LEN * DIM];
    int n = p->len;
    int i, step, energy, new_energy, moved;
    double t;

    //start from a straight chain along the x axis
    for (i = 0; i < n; i++)
    {
        cur[2 * i] = i;
        cur[2 * i + 1] = 0;
    }
    energy = compute_energy(p, cur);

    result.best_energy = energy;
    memcpy(result.coords, cur, 2 * n * sizeof(int));
    result.time_moves = 0.0;
    result.time_energy = 0.0;

    for (step = 0; step < num_steps; step++)
    {
        //propose a move on a copy of the current fold
        t = getTime();
        memcpy(trial, cur, 2 * n * sizeof(int));
        moved = pullMove(cur, trial, n, &seed);
        result.time_moves += getTime() - t;
        if (!moved)
            continue;   //no valid move this step

        t = getTime();
        new_energy = compute_energy(p, trial);
        result.time_energy += getTime() - t;

        //Metropolis rule: always keep a better fold, sometimes keep a worse one
        if (new_energy <= energy ||
            rand_r(&seed) / (RAND_MAX + 1.0) < exp(-(new_energy - energy) / TEMPERATURE))
        {
            memcpy(cur, trial, 2 * n * sizeof(int));
            energy = new_energy;
            if (energy < result.best_energy)
            {
                result.best_energy = energy;
                memcpy(result.coords, cur, 2 * n * sizeof(int));
            }
        }
    }
    return result;
}

//energy = -(number of H-H pairs that touch but are not chain neighbors)
int compute_energy(const Protein *p, const int *coords)
{
    int i, j, contacts = 0;
    for (i = 0; i < p->len; i++)
    {
        if (p->seq[i] != 'H')
            continue;
        for (j = i + 2; j < p->len; j++)   //j = i+1 is a chain neighbor, skip it
        {
            if (p->seq[j] == 'H' &&
                isAdjacent(coords[2 * i], coords[2 * i + 1], coords[2 * j], coords[2 * j + 1]))
                contacts++;
        }
    }
    return -contacts;
}

//read "name sequence best_known_energy" lines; returns how many were read
//(moved here from Srinath's mpi.c so every version shares it)
int parse_proteins(const char *filename, Protein *proteins, int max_proteins)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Cannot open %s\n", filename);
        return 0;
    }
    int count = 0;
    while (count < max_proteins &&
           fscanf(fp, "%31s %127s %d",
                  proteins[count].name,
                  proteins[count].seq,
                  &proteins[count].best_known_energy) == 3)
    {
        proteins[count].len = (int)strlen(proteins[count].seq);
        count++;
    }
    fclose(fp);
    return count;
}

//returns 1 if the fold is a connected chain with no two residues on one square
int check_fold(const Protein *p, const int *coords)
{
    int i, j;
    for (i = 1; i < p->len; i++)
    {
        if (!isAdjacent(coords[2 * i - 2], coords[2 * i - 1], coords[2 * i], coords[2 * i + 1]))
            return 0;   //chain is broken
    }
    for (i = 0; i < p->len; i++)
    {
        for (j = i + 1; j < p->len; j++)
        {
            if (coords[2 * i] == coords[2 * j] && coords[2 * i + 1] == coords[2 * j + 1])
                return 0;   //two residues on the same square
        }
    }
    return 1;
}

//wall clock time in seconds
static double getTime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

//returns 1 if some residue of the chain sits on (x, y)
static int isOccupied(const int *coords, int n, int x, int y)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (coords[2 * i] == x && coords[2 * i + 1] == y)
            return 1;
    }
    return 0;
}

//returns 1 if the two squares share an edge
static int isAdjacent(int x1, int y1, int x2, int y2)
{
    return abs(x1 - x2) + abs(y1 - y2) == 1;
}

/*
 * Try one random pull move. cur is the current fold, trial starts as a
 * copy of it and gets the new fold. Returns 1 if a move was made, 0 if
 * the random choice was blocked.
 *
 * Pick residue i and a direction d (+1 or -1). Residue i+d stays put and
 * residues i-d, i-2d, ... follow behind i:
 *   - i moves to a free square L next to i+d (diagonal to i)
 *   - i-d moves to the square C next to both L and i's old square
 *   - each next residue j takes the old square of residue j+2d, until
 *     one is already next to the residue ahead of it
 * If i is an end of the chain (i+d is off the chain), i moves two free
 * squares away instead and the rest of the chain follows the same way.
 * The new fold is always connected and self-avoiding.
 */
static int pullMove(const int *cur, int *trial, int n, unsigned int *seed)
{
    int i = rand_r(seed) % n;
    int d = (rand_r(seed) % 2) ? 1 : -1;
    int lx, ly, cx, cy, j, k;

    if (i + d < 0 || i + d >= n)
    {
        //end move: C is next to i, L is next to C
        k = rand_r(seed) % 4;
        cx = cur[2 * i] + DX[k];
        cy = cur[2 * i + 1] + DY[k];
        k = rand_r(seed) % 4;
        lx = cx + DX[k];
        ly = cy + DY[k];
        if (isOccupied(cur, n, cx, cy) || isOccupied(cur, n, lx, ly))
            return 0;
    }
    else
    {
        //turn the step from i to i+d by 90 degrees, either way
        int vx = cur[2 * (i + d)] - cur[2 * i];
        int vy = cur[2 * (i + d) + 1] - cur[2 * i + 1];
        int px = -vy, py = vx;
        if (rand_r(seed) % 2)
        {
            px = vy;
            py = -vx;
        }
        lx = cur[2 * (i + d)] + px;
        ly = cur[2 * (i + d) + 1] + py;
        cx = cur[2 * i] + px;
        cy = cur[2 * i + 1] + py;
        if (isOccupied(cur, n, lx, ly))
            return 0;

        //i is a chain end, or i-d already sits on C: only i has to move
        if (i - d < 0 || i - d >= n ||
            (cur[2 * (i - d)] == cx && cur[2 * (i - d) + 1] == cy))
        {
            trial[2 * i] = lx;
            trial[2 * i + 1] = ly;
            return 1;
        }
        if (isOccupied(cur, n, cx, cy))
            return 0;
    }

    //move i to L and i-d to C, then pull the rest of the chain behind them
    trial[2 * i] = lx;
    trial[2 * i + 1] = ly;
    trial[2 * (i - d)] = cx;
    trial[2 * (i - d) + 1] = cy;
    for (j = i - 2 * d; j >= 0 && j < n; j -= d)
    {
        if (isAdjacent(cur[2 * j], cur[2 * j + 1], trial[2 * (j + d)], trial[2 * (j + d) + 1]))
            break;   //chain is connected again, the rest stays put
        trial[2 * j] = cur[2 * (j + 2 * d)];
        trial[2 * j + 1] = cur[2 * (j + 2 * d) + 1];
    }
    return 1;
}
