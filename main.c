#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "cpdb.h"
#include "saxs.h"

void write_saxs(FILE* dest, const double* q, const double* i, const double* e, const size_t len)
{
    if(e != NULL)
    {
        for(size_t x = 0; x < len; x++)
        {
            fprintf(dest, "%.5f\t%f\t%f\n", q[x], i[x], e[x]);
        }
    
    }
    else
    {
        for(size_t x = 0; x < len; x++)
        {
            fprintf(dest, "%.5f\t%f\n", q[x], i[x]);
        }
    }
    fputc('\n', dest);
}
    

int main()
{
    /*
     * Read in pdb file and set up structures for SAXS
     */
    char filename[] = "4fcy.pdb";
    
    size_t buffsize = 500;
    Atom* atoms = malloc(sizeof(Atom) * buffsize);
    size_t atoms_read = read_atoms_from_pdbfile(&atoms, buffsize, fopen(filename,"r"), false);
    printf("Read %lu atoms from %s\n", (unsigned long)atoms_read, filename);
    
    printf("Determining bonding configuration..."); fflush(stdout);
    assign_organic_atom_types(atoms, atoms_read);
    printf("done\n"); fflush(stdout);
    printf("Setting up neighbor lists..."); fflush(stdout);
    NeighborList nl;
    allocate_make_neighbor_list(&nl, atoms, atoms_read, 3.0 + 3.0 + 1.4 * 2);
    printf("done\n"); fflush(stdout);
    printf("Calculating SASA..."); fflush(stdout);
    double* per_atom_sasa = malloc(sizeof(double) * atoms_read);
    shrake_rupley_sasa(per_atom_sasa, 1.4, atoms, atoms_read, true, &nl, NULL, 99);
    printf("done\n"); fflush(stdout);
    
    static const size_t Q_POINTS = 501;
    static const double DELTA_Q = 0.001;
    static const size_t GV_POINTS = 149;
    static const double WATERDEN = 0.334;
    static const size_t RADIUS_SAMPLES = 1;
    static const double RADIUS_MIN = 0.95;
    static const double RADIUS_MAX = 1.05;
    
    printf("Initializing sampling vectors..."); fflush(stdout);
    coord sampling_vectors[GV_POINTS];
    golden_spiral_distribution(sampling_vectors, GV_POINTS);
    printf("done\n"); fflush(stdout);
    
    printf("Initializing form factors..."); fflush(stdout);
    double* q = malloc(sizeof(double) * Q_POINTS);
    for(size_t i = 0; i < Q_POINTS; i++)
        q[i] = DELTA_Q * i;

    double* intensities = malloc(sizeof(double) * Q_POINTS);
    SaxsProfile saxsProfile = { .q = q, .i = intensities, .e = NULL, .q_spacing = DELTA_Q, .length = Q_POINTS};

    FormFactorTable formFactorTable;
    allocate_form_factor_table(&formFactorTable, Q_POINTS, RADIUS_SAMPLES);
    populate_form_factor_table(&formFactorTable, atoms, atoms_read, q, Q_POINTS, RADIUS_MIN, RADIUS_MAX, RADIUS_SAMPLES, WATERDEN);
    printf("done\n"); fflush(stdout);
    
    /*
     * Do the SAXS calculation
     */
    printf("Calculating unfitted scattering..."); fflush(stdout);
    unfitted_saxs(&saxsProfile, atoms, per_atom_sasa, atoms_read, sampling_vectors, GV_POINTS, &formFactorTable, 0, 2.0);
    printf("done\n"); fflush(stdout);
    
    char outfilename[256];
    sprintf(outfilename, "%s%s", filename, "-fastcat.dat");
    printf("Writing unfitted scatting to %s...", outfilename); fflush(stdout);
    FILE* outfile = fopen(outfilename, "w");
    if(outfile != NULL)
    {
        write_saxs(outfile, q, intensities, NULL, Q_POINTS);
        fclose(outfile);
    }
    printf("done\n");
    
    /*
    unfitted_debeye(&saxsProfile, atoms, per_atom_sasa, atoms_read, &formFactorTable, 0, 2.0);
    
    sprintf(outfilename, "%s%s", filename, "-debeye.dat");
    outfile = fopen(outfilename, "w");
    if(outfile != NULL)
    {
        write_saxs(outfile, q, intensities, NULL, Q_POINTS);
        fclose(outfile);
    }
     */

    free_neighborlist(&nl);
    free(per_atom_sasa);
    free(atoms);
    free(q);
    free(intensities);
    free(formFactorTable.vacuo_form_factors);
    free(formFactorTable.dummy_form_factors);
    free(formFactorTable.solvent_form_factors);
    
    return 0;
}