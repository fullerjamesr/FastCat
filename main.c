#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "cpdb.h"
#include "saxs.h"
#include "argparse.h"

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


int main(int argc, const char** argv)
{
    /*
     * Usage
     */
    const char usage[] = "fastcat [options] structure.pdb [experimental scattering.dat]";
    const char epilog[] = "fastcat is developed by James Fuller (fullerj@uchicago.edu)";
    
    /*
     * Arguments
     */
    bool QUIET = false;
    char* mode = "fastcat";
    int Q_POINTS = 501;
    double DELTA_Q = 0.001;
    int GV_POINTS = 149;
    double SOLVENT_DENSITY = 0.334;
    int RADIUS_SAMPLES = 10;
    double RADIUS_MIN = 0.95;
    double RADIUS_MAX = 1.05;
    int HYDRATION_SAMPLES = 30;
    double HYDRATION_MIN = -2.0;
    double HYDRATION_MAX = 4.0;
#ifdef _OPENMP
    int NUM_THREADS = omp_get_num_procs();
#endif
    
    ArgumentOption options[] =
    {
        OPT_HELP("print usage and exit"),
        OPT_BOOLEAN(0, "quiet", &QUIET, "suppress standard output"),
        OPT_STRING('m', "mode", &mode, "method to compute scattering (`fastcat` or `debeye`"),
        OPT_INTEGER('q', "profile_points", &Q_POINTS, "number of points in the calculated scattering profile"),
        OPT_DOUBLE('d', "delta_q", &DELTA_Q, "spacing between points in the calculated scattering profile"),
        OPT_INTEGER('n', "sampling_points", &GV_POINTS, "number of scattering vectors to sample during profile calculation"),
        OPT_DOUBLE('s', "solvent_density", &SOLVENT_DENSITY, "solvent density to use when calculating the scattering contribution from excluded volume"),
        OPT_INTEGER(0, "hydration_layer_samples", &HYDRATION_SAMPLES, "number of hydration layer density samples to enumerate between `hydration_layer_min` and `hydration_layer_max"),
        OPT_DOUBLE(0, "hydration_layer_min", &HYDRATION_MIN, "minimum number of hydration blobs neighboring a fully exposed atom"),
        OPT_DOUBLE(0, "hydration_layer_max", &HYDRATION_MAX, "maximum number of hydration blobs neighboring a fully exposed atom"),
        OPT_INTEGER(0, "excluded_radius_samples", &RADIUS_SAMPLES, "number of adjustments to dummy atom radii to enumerate between `radius_min` and `radius_max`"),
        OPT_DOUBLE(0, "excluded_radius_min", &RADIUS_MIN, "minimum adjustment factor to dummy atom radii"),
        OPT_DOUBLE(0, "excluded_radius_max", &RADIUS_MAX, "maximum adjustment factor to dummy atom radii"),
#ifdef _OPENMP
        OPT_INTEGER('t', "num_threads", &NUM_THREADS, "number of threads to spawn during profile calculation"),
#endif
        OPT_END()
    };
    
    ArgParseInfo argParseInfo = { options, usage, NULL, epilog};
    argc = do_argparse(&argParseInfo, argc, argv);
    if(!QUIET)
    {
        print_state(&argParseInfo);
        putc('\n', stdout);
    }
    
#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#endif
    
    /*
     * Read in pdb file and set up structures for SAXS
     */
    const char* filename = "4fcy.pdb";
    for(size_t i = 0; i < (size_t) argc; i++)
    {
        const char* token = argv[i];
        char* dot_position = strrchr(token, '.');
        if(dot_position && !strcmp(dot_position, ".pdb"))
        {
            filename = token;
        }
        // TODO: looking for input scattering files goes here
        // TODO: if neither a pdb or valid 2/3 column file is found, exit(1) with appropriate error message
    }
    FILE* fh = fopen(filename, "r");
    
    size_t buffsize = 500;
    Atom* atoms = malloc(sizeof(Atom) * buffsize);
    size_t atoms_read = read_atoms_from_pdbfile(&atoms, buffsize, fh, false);
    if(!QUIET)
    {
        printf("Read %lu atoms from %s\n", (unsigned long) atoms_read, filename);
        printf("Determining bonding configuration...");
        fflush(stdout);
    }
    assign_organic_atom_types(atoms, atoms_read);
    if(!QUIET)
    {
        printf("done\n");
        printf("Setting up neighbor lists...");
        fflush(stdout);
    }
    NeighborList nl;
    allocate_make_neighbor_list(&nl, atoms, atoms_read, 3.0 + 3.0 + 1.4 * 2);
    if(!QUIET)
    {
        printf("done\n");
        printf("Calculating SASA...");
        fflush(stdout);
    }
    double* per_atom_sasa = malloc(sizeof(double) * atoms_read);
    shrake_rupley_sasa(per_atom_sasa, 1.4, atoms, atoms_read, true, &nl, NULL, 99);
    if(!QUIET)
    {
        printf("done\n");
        printf("Initializing sampling vectors...");
        fflush(stdout);
    }
    coord sampling_vectors[GV_POINTS];
    golden_spiral_distribution(sampling_vectors, (size_t)GV_POINTS);
    if(!QUIET)
    {
        printf("done\n");
        printf("Initializing form factors...");
        fflush(stdout);
    }
    double* q = malloc(sizeof(double) * Q_POINTS);
    for(size_t i = 0; i < Q_POINTS; i++)
        q[i] = DELTA_Q * i;

    double* intensities = malloc(sizeof(double) * Q_POINTS);
    SaxsProfile saxsProfile = { .q = q, .i = intensities, .e = NULL, .q_spacing = DELTA_Q, .length = (size_t)Q_POINTS};

    FormFactorTable formFactorTable;
    allocate_form_factor_table(&formFactorTable, (size_t)Q_POINTS, (size_t)RADIUS_SAMPLES);
    populate_form_factor_table(&formFactorTable, atoms, atoms_read, q, (size_t)Q_POINTS, RADIUS_MIN, RADIUS_MAX, (size_t)RADIUS_SAMPLES, SOLVENT_DENSITY);
    if(!QUIET)
    {
        printf("done\n");
        printf("Calculating unfitted scattering...");
        fflush(stdout);
    }
    unfitted_saxs(&saxsProfile, atoms, per_atom_sasa, atoms_read, sampling_vectors, (size_t)GV_POINTS, &formFactorTable, 0, 2.0);
    char outfilename[256];
    sprintf(outfilename, "%s%s", filename, "-fastcat.dat");
    if(!QUIET)
    {
        printf("done\n");
        printf("Writing unfitted scatting to %s...", outfilename);
        fflush(stdout);
    }
    FILE* outfile = fopen(outfilename, "w");
    if(outfile != NULL)
    {
        write_saxs(outfile, q, intensities, NULL, (size_t)Q_POINTS);
        fclose(outfile);
    }
    if(!QUIET)
    {
        printf("done\n");
    }
    
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