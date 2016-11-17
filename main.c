#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include "cpdb.h"
#include "saxs.h"

int main()
{
    char filename[] = "4fcy.pdb";
    
    size_t buffsize = 500;
    Atom* atoms = malloc(sizeof(Atom) * buffsize);
    size_t atoms_read = read_atoms_from_pdbfile(&atoms,buffsize,fopen(filename,"r"),false);
    assign_organic_atom_types(atoms,atoms_read);
    
    
    static const size_t Q_POINTS = 201;
    static const double DELTA_Q = 0.01;
    static const size_t GV_POINTS = 149;
    static const double waterden = 0.334;
    
    coord sampling_vectors[GV_POINTS];
    golden_spiral_distribution(sampling_vectors,GV_POINTS);
    
    double q[Q_POINTS];
    for(size_t i = 0; i < Q_POINTS; i++)
        q[i] = DELTA_Q * i;
    
    double intensities[Q_POINTS];
    
    clock_t begin, end;
    double time_spent;
    begin = clock();
    unfitted_saxs(intensities, atoms, atoms_read, q, Q_POINTS, sampling_vectors, GV_POINTS, waterden, 1.0);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Original method took %0.4f sec\n", time_spent);
    
    /*
    begin = clock();
    unfitted_saxs2(intensities, atoms, atoms_read, q, Q_POINTS, sampling_vectors, GV_POINTS, waterden, 1.0);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Modified method took %0.4f sec\n", time_spent);
     */
    
    free(atoms);
    return 0;
}