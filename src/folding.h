#ifndef FOLDING_H
#define FOLDING_H

#define DIM 2
#define MAX_SEQ_LEN 128

typedef struct {
    char seq[MAX_SEQ_LEN];
    int  len;
    char name[32];
    int  best_known_energy;
} Protein;

typedef struct {
    int best_energy;
    int coords[MAX_SEQ_LEN * DIM];
} FoldResult;

FoldResult fold_protein(const Protein *p, int num_steps, unsigned int seed);
int        compute_energy(const Protein *p, const int *coords);
int        parse_proteins(const char *filename, Protein *proteins, int max_proteins);

#endif
