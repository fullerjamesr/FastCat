#define __USE_MINGW_ANSI_STDIO 1
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

#ifndef EXPERIMENTAL_DATA_ALLOC_SIZE
    #define EXPERIMENTAL_DATA_ALLOC_SIZE 1024
#endif
#ifndef ATOMS_ALLOC_SIZE
    #define ATOMS_ALLOC_SIZE 10000
#endif

void write_saxs(FILE* destination, const SaxsProfile* const source)
{
    if(source->e)
        write_n_column_data(destination, source->length, 3, source->q, source->i, source->e);
    else
        write_n_column_data(destination, source->length, 2, source->q, source->i);
}


int main(int argc, const char** argv)
{
    /*
     *
     * Usage
     *
     */
    const char usage[] = "fastcat [options] structure.pdb [experimental scattering.dat]";
    const char epilog[] = "fastcat is developed by James Fuller https://github.com/fullerjamesr";
    
    /*
     *
     * Arguments
     *
     */
    bool QUIET = false;
    char* MODE = "fastcat";
    bool OFFSET = false;
    int Q_POINTS = 501;
    double DELTA_Q = 0.001;
    int GV_POINTS = 149;
    double SOLVENT_DENSITY = 0.334;
    int RADIUS_SAMPLES = 20;
    double RADIUS_MIN = 0.95;
    double RADIUS_MAX = 1.05;
    int HYDRATION_SAMPLES = 60;
    double HYDRATION_MIN = -2.0;
    double HYDRATION_MAX = 4.0;
    int SASA_POINTS = 99;
#ifdef _OPENMP
    char* threads_str = getenv("OMP_NUM_THREADS");
    int NUM_THREADS = threads_str ? abs(atoi(threads_str)) : omp_get_num_procs();
#endif
    
    ArgumentOption options[] =
    {
        OPT_HELP("print usage and exit"),
        OPT_BOOLEAN(0, "quiet", &QUIET, "suppress standard output"),
        OPT_STRING('m', "mode", &MODE, "method to compute scattering (one of `fastcat` or `fastcat-search` or `debeye`"),
        OPT_BOOLEAN('o', "offset", &OFFSET, "when fitting profiles, fit a constant offset added to the experimental data to improve the fit"),
        OPT_INTEGER('q', "profile_points", &Q_POINTS, "number of points in the calculated scattering profile"),
        OPT_DOUBLE('d', "delta_q", &DELTA_Q, "spacing between points in the calculated scattering profile"),
        OPT_INTEGER('n', "sampling_vectors", &GV_POINTS, "number of scattering vectors to sample during profile calculation"),
        OPT_DOUBLE('s', "solvent_density", &SOLVENT_DENSITY, "solvent density to use when calculating the scattering contribution from excluded volume"),
        OPT_INTEGER(0, "hydration_layer_samples", &HYDRATION_SAMPLES, "number of hydration layer density samples to enumerate between `hydration_layer_min` and `hydration_layer_max`"),
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
    
    // Make sure radius samples is <= 64
    if(RADIUS_SAMPLES > 64)
        RADIUS_SAMPLES = 64;
    
    /*
     *
     * Read in PDB and (potential) data file.
     * The data file, if present, overrides any command line arguments for the q-grid
     *
     */
    // In theory, if I were to ever support multiple PDB or data files at once, sorting out those files would go below.
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
    if(!pdb_file_handle)
    {
        perror("Could not open PDB file:");
        exit(1);
    }
    FILE* dat_file_handle = fopen(dat_filename, "r");
    
    double* experimental_q = NULL;
    double* experimental_i = NULL;
    double* experimental_e = NULL;
    size_t real_profile_size = (size_t)Q_POINTS;
    if(dat_file_handle)
    {
        size_t profile_alloc_size = EXPERIMENTAL_DATA_ALLOC_SIZE;
        experimental_q = malloc(sizeof(double) * profile_alloc_size);
        experimental_i = malloc(sizeof(double) * profile_alloc_size);
        experimental_e = malloc(sizeof(double) * profile_alloc_size);
        
        double* data[3] = { experimental_q, experimental_i, experimental_e};
        real_profile_size = read_3_column_data(data, profile_alloc_size, dat_file_handle);
        
        // read until there are no more lines coming in
        while(real_profile_size == profile_alloc_size)
        {
            profile_alloc_size *= 2;
            double* re_q = realloc(experimental_q, sizeof(double) * profile_alloc_size);
            double* re_i = realloc(experimental_i, sizeof(double) * profile_alloc_size);
            double* re_e = realloc(experimental_e, sizeof(double) * profile_alloc_size);
            if(!re_q || !re_i || !re_e)
                break;
            
            experimental_q = re_q;
            experimental_i = re_i;
            experimental_e = re_e;
            data[0] = experimental_q + real_profile_size;
            data[1] = experimental_i + real_profile_size;
            data[2] = experimental_e + real_profile_size;
            real_profile_size += read_3_column_data(data, profile_alloc_size - real_profile_size, dat_file_handle);
        }
        
        fclose(dat_file_handle);
        
        if(!QUIET)
            printf("Read %lu experimental data points from %s\n", (unsigned long)real_profile_size, dat_filename);
       
        Q_POINTS = (int)real_profile_size;
        DELTA_Q = experimental_q[1] - experimental_q[0];
    }
    else
    {
        experimental_q = malloc(sizeof(double) * Q_POINTS);
        experimental_i = malloc(sizeof(double) * Q_POINTS);
        for(size_t i = 0; i < Q_POINTS; i++)
            experimental_q[i] = DELTA_Q * i;
    }
    SaxsProfile experimental_data = { .q = experimental_q, .i = experimental_i, .e = experimental_e, .length = real_profile_size };
    
    // Printing the status of all options can finally be done here because we know if the Q_GRID is specified by experimental data
    print_state(&argParseInfo);
    
    size_t atoms_buffersize = ATOMS_ALLOC_SIZE;
    Atom* atoms = malloc(sizeof(Atom) * atoms_buffersize);
    size_t atoms_read = read_atoms_from_pdbfile(&atoms, atoms_buffersize, pdb_file_handle, false);
    fclose(pdb_file_handle);
    if(!QUIET)
    {
        printf("\nRead %lu atoms from %s\n", (unsigned long) atoms_read, pdb_filename);
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
    allocate_make_neighbor_list(&nl, atoms, atoms_read, 9.0);
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
    
    FormFactorTable formFactorTable;
    allocate_form_factor_table(&formFactorTable, real_profile_size, (size_t)RADIUS_SAMPLES);
    populate_form_factor_table(&formFactorTable, atoms, atoms_read, experimental_q, real_profile_size, RADIUS_MIN, RADIUS_MAX, (size_t)RADIUS_SAMPLES, SOLVENT_DENSITY);
    if(!QUIET)
        puts("done");
    
    SaxsProfile calculated_profile;
    if(dat_file_handle)
    {
        if(!QUIET)
        {
            fputs("Enumerating scattering profiles...", stdout);
            fflush(stdout);
        }
        double* intensities = malloc(sizeof(double) * experimental_data.length * RADIUS_SAMPLES * HYDRATION_SAMPLES);
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
        fitted_saxs(intensities, experimental_q, experimental_data.length, atoms, per_atom_sasa, atoms_read, sampling_vectors, (size_t)GV_POINTS, &formFactorTable, HYDRATION_MIN, HYDRATION_MAX, (size_t)HYDRATION_SAMPLES);

        size_t best_d = 0;
        size_t best_h = 0;
        double best_chisq = DBL_MAX;
        ProfileScaleResult best_scale = { .constant = 0.0, .offset = 0.0 };
        for(size_t d = 0; d < RADIUS_SAMPLES; d++)
        {
            for(size_t h = 0; h < HYDRATION_SAMPLES; h++)
            {
                SaxsProfile* model = profiles + d * HYDRATION_SAMPLES + h;
                ProfileScaleResult scale_factors = OFFSET
                                                   ? find_chisq_scale_factor_with_offset(&experimental_data, model)
                                                   : find_chisq_scale_factor(&experimental_data, model);
                for(size_t i = 0; i < experimental_data.length; i++)
                    model->i[i] = model->i[i] * scale_factors.constant - scale_factors.offset;
                double score = chi_square(&experimental_data, model);
                if(score < best_chisq)
                {
                    best_d = d;
                    best_h = h;
                    best_chisq = score;
                    best_scale = scale_factors;
                }
            }
        }
        
        if(!QUIET)
            printf("done\n\nBest fit:\ndummy atom radius adjustment = %f\nhydration density factor = %f\n"
                   "scaling factor = %e\nprofile offset = %f\nchisquare = %f\n\n",
                   RADIUS_MIN + (RADIUS_MAX - RADIUS_MIN) * best_d / RADIUS_SAMPLES,
                   HYDRATION_MIN + (HYDRATION_MAX - HYDRATION_MIN) * best_h / HYDRATION_SAMPLES,
                   best_scale.constant, -best_scale.offset, best_chisq);
        calculated_profile.q = experimental_data.q;
        calculated_profile.i = profiles[best_d * HYDRATION_SAMPLES + best_h].i;
        calculated_profile.e = NULL;
        calculated_profile.length = experimental_data.length;
        
        free(profiles);
    }
    else
    {
        if(!QUIET)
        {
            fputs("Calculating scattering...", stdout);
            fflush(stdout);
        }
        calculated_profile.q = experimental_q;
        calculated_profile.i = experimental_i;
        calculated_profile.e = NULL;
        calculated_profile.length = (size_t)Q_POINTS;
        unfitted_saxs(&calculated_profile, atoms, per_atom_sasa, atoms_read, sampling_vectors, (size_t) GV_POINTS,
                      &formFactorTable, (size_t)RADIUS_SAMPLES/2, HYDRATION_MIN + (HYDRATION_MAX - HYDRATION_MIN) / 2);
        puts("done");
    }
    
    
    char outfilename[256];
    sprintf(outfilename, "%s%s", pdb_filename, "-fastcat.dat");
    FILE* outfile = fopen(outfilename, "w");
    if(outfile)
    {
        printf("Writing result to %s...", outfilename);
        write_saxs(outfile, &calculated_profile);
        fclose(outfile);
        puts("done");
    }
    else
        perror("Could not write resulting profile:");
    
    
    free_neighborlist(&nl);
    free(per_atom_sasa);
    free(atoms);
    free(experimental_i);
    if(experimental_e)
        free(experimental_e);
    free(formFactorTable.vacuo_form_factors);
    free(formFactorTable.dummy_form_factors);
    free(formFactorTable.solvent_form_factors);
    
    return 0;
}