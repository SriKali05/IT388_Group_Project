# IT388_Group_Project

## Input file: data/proteins.txt

One protein per line: name sequence best_known_energy

    S20 HPHPPHHPHPPHPHHPPHPH -9

- name: label (S20 = 20 residues)
- sequence: H = hydrophobic, P = polar
- best_known_energy: lowest known energy on the 2D square lattice. Used only
  to check our results, not for folding.

No header, comments, or blank lines. Sequences are the standard 2D HP
benchmarks (S1-1 to S1-8) from Table 1 of ref [1], expanded by hand from the
paper's compressed notation. Energies apply to 2D only (DIM 2).

Groups for the problem-size study:
    Short:  S20, S24, S25
    Medium: S36, S48, S50
    Long:   S60, S64