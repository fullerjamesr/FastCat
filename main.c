#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <float.h>
#include "cpdb.h"
#include "saxs.h"
#include "argparse.h"
#include "fileio.h"

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
    char* MODE = "fastcat";
    int Q_POINTS = 501;
    double DELTA_Q = 0.001;
    int GV_POINTS = 149;
    double SOLVENT_DENSITY = 0.334;
    int RADIUS_SAMPLES = 10;
    double RADIUS_MIN = 0.95;
    double RADIUS_MAX = 1.05;
    int HYDRATION_SAMPLES = 15;
    double HYDRATION_MIN = -2.0;
    double HYDRATION_MAX = 4.0;
    int SASA_POINTS = 99;
#ifdef _OPENMP
    int NUM_THREADS = omp_get_num_procs();
#endif
    
    ArgumentOption options[] =
    {
        OPT_HELP("print usage and exit"),
        OPT_BOOLEAN(0, "quiet", &QUIET, "suppress standard output"),
        OPT_STRING('m', "mode", &MODE, "method to compute scattering (one of `fastcat` or `fastcat-search` or `debeye`"),
        OPT_INTEGER('q', "profile_points", &Q_POINTS, "number of points in the calculated scattering profile"),
        OPT_DOUBLE('d', "delta_q", &DELTA_Q, "spacing between points in the calculated scattering profile"),
        OPT_INTEGER('n', "sampling_vectors", &GV_POINTS, "number of scattering vectors to sample during profile calculation"),
        OPT_DOUBLE('s', "solvent_density", &SOLVENT_DENSITY, "solvent density to use when calculating the scattering contribution from excluded volume"),
        OPT_INTEGER(0, "hydration_layer_samples", &HYDRATION_SAMPLES, "number of hydration layer density samples to enumerate between `hydration_layer_min` and `hydration_layer_max"),
        OPT_DOUBLE(0, "hydration_layer_min", &HYDRATION_MIN, "minimum number of hydration blobs neighboring a fully exposed atom"),
        OPT_DOUBLE(0, "hydration_layer_max", &HYDRATION_MAX, "maximum number of hydration blobs neighboring a fully exposed atom"),
        OPT_INTEGER(0, "excluded_radius_samples", &RADIUS_SAMPLES, "number of adjustments to dummy atom radii to enumerate between `radius_min` and `radius_max`"),
        OPT_DOUBLE(0, "excluded_radius_min", &RADIUS_MIN, "minimum adjustment factor to dummy atom radii"),
        OPT_DOUBLE(0, "excluded_radius_max", &RADIUS_MAX, "maximum adjustment factor to dummy atom radii"),
        OPT_INTEGER(0, "sasa_sampling_vectors", &SASA_POINTS, "number of test vectors to use to calculate SASA for hydration layer"),
#ifdef _OPENMP
        OPT_INTEGER('t', "num_threads", &NUM_THREADS, "number of threads to spawn during profile calculation"),
#endif
        OPT_END()
    };
    
    ArgParseInfo argParseInfo = { options, usage, NULL, epilog};
    argc = do_argparse(&argParseInfo, argc, argv);
    
#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#endif
    
    /*
     * Read in PDB file and [potential] .dat file
     */
    const char* pdb_filename = "4fcy.pdb";
    const char* dat_filename = "4fcy.dat";
    for(size_t i = 0; i < (size_t) argc; i++)
    {
        const char* token = argv[i];
        char* dot_position = strrchr(token, '.');
        if(dot_position && !strcmp(dot_position, ".pdb"))
        {
            pdb_filename = token;
        }
        else
        {
            dat_filename = token;
        }
    }
    
    FILE* pdb_file_handle = fopen(pdb_filename, "r");
    FILE* dat_file_handle = fopen(dat_filename, "r");
    if(!pdb_file_handle)
    {
        fprintf(stderr, "error: could not open file %s\n", pdb_filename);
        exit(1);
    }
    
    double* experimental_q = NULL;
    double* experimental_i = NULL;
    double* experimental_e = NULL;
    SaxsProfile experimental_data;
    if(dat_file_handle)
    {
        experimental_q = malloc(sizeof(double) * 1024);
        experimental_i = malloc(sizeof(double) * 1024);
        experimental_e = malloc(sizeof(double) * 1024);
        experimental_data.q = experimental_q;
        experimental_data.i = experimental_i;
        experimental_data.e = experimental_e;
        experimental_data.length = 1024;
        double* data[3] = { experimental_q, experimental_i, experimental_e};
        size_t real_size = read_3_column_data(data, 1024, dat_file_handle);
        experimental_data.length = real_size;
        Q_POINTS = (int) real_size;
        fclose(dat_file_handle);
        if(!QUIET)
            printf("Read %lu experimental data points from %s\n", (unsigned long)real_size, dat_filename);
    }
    
    size_t buffsize = 500;
    Atom* atoms = malloc(sizeof(Atom) * buffsize);
    size_t atoms_read = read_atoms_from_pdbfile(&atoms, buffsize, pdb_file_handle, false);
    fclose(pdb_file_handle);
    if(!QUIET)
    {
        printf("Read %lu atoms from %s\n", (unsigned long) atoms_read, pdb_filename);
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
    shrake_rupley_sasa(per_atom_sasa, 1.4, atoms, atoms_read, true, &nl, NULL,(size_t)SASA_POINTS);
    if(!QUIET)
    {
        puts("done");
        fputs("Initializing sampling vectors...", stdout);
        fflush(stdout);
    }
    coord sampling_vectors[GV_POINTS];
    golden_spiral_distribution(sampling_vectors, (size_t)GV_POINTS);
    if(!QUIET)
    {
        puts("done");
        fputs("Initializing form factors...", stdout);
        fflush(stdout);
    }
    
    double* q;
    if(dat_file_handle)
        q = experimental_q;
    else
    {
        q = malloc(sizeof(double) * Q_POINTS);
        for(size_t i = 0; i < Q_POINTS; i++)
            q[i] = DELTA_Q * i;
    }
    FormFactorTable formFactorTable;
    allocate_form_factor_table(&formFactorTable, (size_t)Q_POINTS, (size_t)RADIUS_SAMPLES);
    populate_form_factor_table(&formFactorTable, atoms, atoms_read, q, (size_t)Q_POINTS, RADIUS_MIN, RADIUS_MAX, (size_t)RADIUS_SAMPLES, SOLVENT_DENSITY);
    if(!QUIET)
        puts("done");
    
    SaxsProfile calculated_profile;
    double* intensities;
    if(dat_file_handle)
    {
        if(!QUIET)
            fputs("Enumerating scattering profiles...", stdout);
        intensities = malloc(sizeof(double) * experimental_data.length * RADIUS_SAMPLES * HYDRATION_SAMPLES);
        double* intensities_marker = intensities;
        SaxsProfile* profiles = malloc(sizeof(SaxsProfile) * RADIUS_SAMPLES * HYDRATION_SAMPLES);
        for(size_t d = 0; d < RADIUS_SAMPLES; d++)
        {
            for(size_t h = 0; h < HYDRATION_SAMPLES; h++)
            {
                SaxsProfile* profile = profiles + d * HYDRATION_SAMPLES + h;
                profile->q = experimental_q;
                profile->i = intensities_marker;
                profile->e = NULL;
                intensities_marker += experimental_data.length;
            }
        }
        fitted_saxs(intensities, q, experimental_data.length, atoms, per_atom_sasa, atoms_read, sampling_vectors, (size_t)GV_POINTS, &formFactorTable, HYDRATION_MIN, HYDRATION_MAX, (size_t)HYDRATION_SAMPLES);

        size_t best_d = 0;
        size_t best_h = 0;
        double best_chisq = DBL_MAX;
        for(size_t d = 0; d < RADIUS_SAMPLES; d++)
        {
            for(size_t h = 0; h < HYDRATION_SAMPLES; h++)
            {
                SaxsProfile* model = profiles + d * HYDRATION_SAMPLES + h;
                ProfileScaleResult scale_factors = find_chisq_scale_factor(&experimental_data, model);
                for(size_t i = 0; i < experimental_data.length; i++)
                    model->i[i] *= scale_factors.constant;
                double score = chi_square(&experimental_data, profiles + d * HYDRATION_SAMPLES + h);
                if(score < best_chisq)
                {
                    best_d = d;
                    best_h = h;
                    best_chisq = score;
                }
            }
        }
        
        if(!QUIET)
            printf("done\nBest chisquare is %f with radius adjustment factor = %f and hydration density = %f\n", best_chisq,
                   RADIUS_MIN + (RADIUS_MAX - RADIUS_MIN) * best_d / RADIUS_SAMPLES,
                   HYDRATION_MIN + (HYDRATION_MAX - HYDRATION_MIN) * best_h / HYDRATION_SAMPLES);
        
        calculated_profile.q = experimental_data.q;
        calculated_profile.i = profiles[best_d * HYDRATION_SAMPLES + best_h].i;
        calculated_profile.e = NULL;
        calculated_profile.length = experimental_data.length;
        
        free(profiles);
    }
    else
    {
        if(!QUIET)
            fputs("Calculating scattering...", stdout);
        intensities = malloc(sizeof(double) * Q_POINTS);
        calculated_profile.q = q;
        calculated_profile.i = intensities;
        calculated_profile.e = NULL;
        calculated_profile.length = (size_t)Q_POINTS;
        unfitted_saxs(&calculated_profile, atoms, per_atom_sasa, atoms_read, sampling_vectors, (size_t) GV_POINTS,
                      &formFactorTable, (size_t)RADIUS_SAMPLES/2, HYDRATION_MIN + (HYDRATION_MAX - HYDRATION_MIN) / 2);
        puts("done");
    }
    
    
    char outfilename[256];
    sprintf(outfilename, "%s%s", pdb_filename, "-fastcat.dat");
    FILE* outfile = fopen(outfilename, "w");
    if(outfile != NULL)
    {
        printf("Writing result to %s...", outfilename);
        write_saxs(outfile, calculated_profile.q, calculated_profile.i, calculated_profile.e, calculated_profile.length);
        fclose(outfile);
        puts("done");
    }
    
    
    free_neighborlist(&nl);
    free(per_atom_sasa);
    free(atoms);
    free(q);
    free(intensities);
    if(experimental_i)
        free(experimental_i);
    if(experimental_e)
        free(experimental_e);
    free(formFactorTable.vacuo_form_factors);
    free(formFactorTable.dummy_form_factors);
    free(formFactorTable.solvent_form_factors);
    
    return 0;
}