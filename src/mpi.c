#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "folding.h"

#define MAX_PROTEINS 64
#define NUM_SEARCHES 1000
#define MC_STEPS     5000

int parse_proteins(const char *filename, Protein *proteins, int max_proteins) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return 0;
    }
    int count = 0;
    while (count < max_proteins &&
           fscanf(fp, "%31s %127s %d",
                  proteins[count].name,
                  proteins[count].seq,
                  &proteins[count].best_known_energy) == 3) {
        proteins[count].len = (int)strlen(proteins[count].seq);
        count++;
    }
    fclose(fp);
    return count;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const char *data_file = "data/proteins.txt";
    if (argc > 1)
        data_file = argv[1];

    Protein proteins[MAX_PROTEINS];
    int num_proteins = 0;

    if (rank == 0) {
        num_proteins = parse_proteins(data_file, proteins, MAX_PROTEINS);
        if (num_proteins == 0) {
            fprintf(stderr, "No proteins loaded.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Loaded %d proteins, running on %d ranks\n", num_proteins, nprocs);
    }

    MPI_Bcast(&num_proteins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(proteins, num_proteins * (int)sizeof(Protein), MPI_BYTE, 0, MPI_COMM_WORLD);

    int my_start = (rank * num_proteins) / nprocs;
    int my_end   = ((rank + 1) * num_proteins) / nprocs;
    int my_count = my_end - my_start;

    int local_energies[MAX_PROTEINS];
    memset(local_energies, 0, sizeof(local_energies));

    double t_start = MPI_Wtime();

    for (int i = 0; i < my_count; i++) {
        int pi = my_start + i;
        int best_energy = 0;

        for (int s = 0; s < NUM_SEARCHES; s++) {
            unsigned int seed = (unsigned int)(rank * 1000 + pi * 100 + s);
            FoldResult result = fold_protein(&proteins[pi], MC_STEPS, seed);
            if (result.best_energy < best_energy)
                best_energy = result.best_energy;
        }

        local_energies[i] = best_energy;
    }

    double t_end = MPI_Wtime();
    double local_time = t_end - t_start;

    int    recv_counts[MAX_PROTEINS];
    int    displs[MAX_PROTEINS];
    int    all_energies[MAX_PROTEINS];

    for (int r = 0; r < nprocs; r++) {
        int rs = (r * num_proteins) / nprocs;
        int re = ((r + 1) * num_proteins) / nprocs;
        recv_counts[r] = re - rs;
        displs[r] = rs;
    }

    MPI_Gatherv(local_energies, my_count, MPI_INT,
                all_energies, recv_counts, displs, MPI_INT,
                0, MPI_COMM_WORLD);

    double max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\n%-6s %6s %10s %10s\n", "Name", "Len", "Best_E", "Known_E");
        printf("--------------------------------------\n");
        for (int i = 0; i < num_proteins; i++) {
            printf("%-6s %6d %10d %10d\n",
                   proteins[i].name,
                   proteins[i].len,
                   all_energies[i],
                   proteins[i].best_known_energy);
        }
        printf("\nWall time: %.4f seconds (%d ranks)\n", max_time, nprocs);
    }

    MPI_Finalize();
    return 0;
}
