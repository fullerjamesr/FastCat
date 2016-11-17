#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "saxs.h"

static const double bohr = 0.52978;

double fast_form_factor(unsigned int z, double q)
{
    return sqrt(exp(-0.23 * q * q)) * (z+1);
}

static const double FormFactorConstants[30][7] =
{
    {2.0,1.0,1.0,1.00E+05,1.0,1.00E+05,1.0},
    {1.994,1.07342,0.64272,0.12568,0.0845,0.14111,0.17253},
    {2.934,0.38714,0.09257,0.07786,0.554,0.01882,0.00128},
    {2.998,0.28152,0.01825,0.04686,0.1465,0.01342,5.80E-04},
    {2.9952,0.2303,0.00808,0.03316,0.06074,0.01083,3.69E-04},
    {2.87,0.18732,0.00503,0.02471,0.03378,0.01421,5.79E-04},
    {3.0808,0.13286,0.2496,0.02926,0.00127,0.00181,7.60E-04},
    {2.5968,0.0939,0.13992,0.02472,0.00135,0.00813,0.01953},
    {2.676,0.08332,0.09631,0.02226,0.00101,0.00534,0.00698},
    {2.6768,0.07474,0.0664,0.02016,8.15E-04,0.004,0.00361},
    {2.4892,0.0539,0.0382,0.01797,7.56E-04,0.00855,0.11457},
    {2.573,0.05214,0.02571,0.01662,5.95E-04,0.00493,0.03657},
    {2.5948,0.04874,0.01839,0.01539,4.96E-04,0.00361,0.02043},
    {2.7396,0.04866,0.01224,0.01441,3.82E-04,0.00192,0.00521},
    {2.8128,0.04704,0.00832,0.01354,3.14E-04,0.00121,0.00187},
    {2.9188,0.0459,0.00561,0.01276,2.54E-04,6.60E-04,5.16E-04},
    {3.0348,0.04484,0.00368,0.0121,2.06E-04,3.10E-04,1.06E-04},
    {3.1704,0.0439,0.00232,0.01151,1.66E-04,1.07E-04,1.19E-05},
    {3.1948,0.0416,0.00176,0.01105,1.48E-04,7.92E-05,6.49E-06},
    {3.1412,0.03902,0.00155,0.01046,1.38E-04,1.06E-04,1.21E-05},
    {2.5908,0.03036,0.00331,0.00964,1.83E-04,0.00112,0.0016},
    {2.5368,0.02814,0.00307,0.00923,1.73E-04,0.00116,0.00172},
    {2.7128,0.02944,0.00212,0.00877,1.36E-04,6.29E-04,4.76E-04},
    {2.7208,0.02866,0.00167,0.00854,1.26E-04,5.43E-04,3.24E-04},
    {2.536,0.02488,0.00207,0.00813,1.32E-04,8.91E-04,9.25E-04},
    {2.58,0.02458,0.00171,0.00785,1.18E-04,7.26E-04,5.92E-04},
    {2.632,0.02438,0.00137,0.00763,1.05E-04,5.73E-04,3.55E-04},
    {2.648,0.02378,0.00123,0.00727,9.46E-05,5.05E-04,2.69E-04},
    {2.648,0.02328,0.00102,0.00712,8.91E-05,4.61E-04,2.12E-04},
    {2.78,0.02398,7.35E-04,0.00677,7.09E-05,2.30E-04,5.28E-05}
};
double form_factor(unsigned int z, double q)
{
    if(z < 31)
    {
        const double* constants = FormFactorConstants[z];
        z += 1;
        double qbohr = bohr * q / 2;    // This already takes into account the 4pi difference between q standards
        double qbohr2 = qbohr * qbohr;
        double a1z_r = pow(constants[1] * z, constants[0]);
        double a2z_r = pow(constants[3] * z, constants[0]);
        double a3z = constants[5] * z;
        double a3z_r = pow(a3z, constants[0]);
        
        double part1 = a1z_r / pow(a1z_r + constants[2] * pow(qbohr, constants[0]), constants[0]);
        
        double part2denominator = pow(2.0, constants[0]) * a2z_r + constants[4] * qbohr2;
        double part2 = a2z_r / part2denominator / part2denominator;
        
        double part3denominator = a3z * a3z * 4 + constants[6] * qbohr2;
        double part3 = a3z_r / part3denominator / part3denominator;
        
        return part1 + part2 + part3;
    }
    else
        return fast_form_factor(z,q);
}

void unfitted_saxs2(double* restrict intensities, Atom* atoms, size_t num_atoms,
                                      double* restrict q_grid, size_t len_grid, coord* sampling_vectors, size_t num_vectors,
                                      double solvent_density, double radius_adjustment)
{
    // Iterate thru atoms to generate a "transposed" array of positions
    coord* positions = malloc(sizeof(coord) * num_atoms);
    //double** newff = malloc(sizeof(double) * len_grid * )
    for(int i = 0; i < num_atoms; i++)
        positions[i] = atoms[i].position;
    
    
    // The form factors for each atom in a vacuum as depend only on q and atom.element
    double* vacuo_form_factors[AtomTypeCount] = { NULL };
    double* form_factors[AtomTypeCount + OrganicAtomConfigCount] = { NULL };
    size_t* offsets = malloc(sizeof(size_t) * num_atoms);   // maps indices of atom --> form_factors
    // Avoid repeated squaring of q when adjusting dummy atom form factors
    double* q_squared = malloc(sizeof(double) * len_grid);
    
    // These properties can be moved out of the q loop for calculating the corrected form factors
    double vol;
    double pivol23;
    double dummy_ff;
    
    clock_t begin, end;
    double time_spent;
    begin = clock();
    // do hydrogens first to add implicitly to combined organic atoms
    //      also populate the q squared lookup table while we're at it
    vacuo_form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    vol = sphere_volume(VdwRadii[Hydrogen] * radius_adjustment);
    pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
    dummy_ff = solvent_density * vol;
    for(size_t q = 0; q < len_grid; q++)
    {
        double q2 = q_grid[q] * q_grid[q];
        double vff = form_factor(Hydrogen, q_grid[q]);
        double ff = vff - dummy_ff * exp(q2 * pivol23);
        vacuo_form_factors[Hydrogen][q] = vff;
        form_factors[Hydrogen][q] = ff;
        q_squared[q] = q2;
    }
    
    for(size_t a = 0; a < num_atoms; a++)
    {
        offsets[a] = atoms[a].organic_type < UnknownBonds ? AtomTypeCount+atoms[a].organic_type : atoms[a].element;
        assert(offsets[a] >= AtomTypeCount);
        // If this is the first occurance of this element
        // ... then we can calculate the vacuo form factors across all q, while also doing the form factors for this atom
        if(vacuo_form_factors[atoms[a].element] == NULL)
        {
            
            vacuo_form_factors[atoms[a].element] = malloc(sizeof(double) * len_grid);
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            //vol = get_atom_volume(atoms+a,true);
            vol = get_atom_volume(atoms + a, true) * radius_adjustment * radius_adjustment * radius_adjustment;
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = form_factor(atoms[a].element, q_grid[q]);
                vacuo_form_factors[atoms[a].element][q] = vff;
                double ff = vff + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q] - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
            //  Else if we've never seen this combined organic atom before
            //  ... just do the corrected form factor calculation using previously stored vacuo factors
        else if(form_factors[offsets[a]] == NULL)
        {
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            vol = get_atom_volume(atoms + a, true) * radius_adjustment * radius_adjustment * radius_adjustment;
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = vacuo_form_factors[atoms[a].element][q] + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q];
                double ff = vff - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
    }
    // Step 1: Place solvent globs
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Populating form factor table took %0.4f sec\n", time_spent);
    
    begin = clock();
    // Step 2: calculate scattering
    #pragma omp parallel for
    for(size_t q = 0; q < len_grid; q++)
    {
        double intensity = 0.0;
        for(size_t v = 0; v < num_vectors; v++)
        {
            double cos_bucket = 0.0;
            double sin_bucket = 0.0;
            coord q_vector = vec_multiply(sampling_vectors[v], q_grid[q]);
            
            for(size_t a = 0; a < num_atoms; a++)
            {
                double corrected_f = form_factors[offsets[a]][q];
                double scat = dot_product(q_vector,positions[a]);
                cos_bucket += corrected_f * cos(scat);
                sin_bucket += corrected_f * sin(scat);
            }
            intensity += cos_bucket * cos_bucket + sin_bucket * sin_bucket;
        }
        intensities[q] = intensity / num_vectors;
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Calculating scattering took %0.4f sec\n", time_spent);
    
    // Free allocated memory
    free(offsets);
    free(q_squared);
    free(positions);
    for(size_t e = 0; e < AtomTypeCount + OrganicAtomConfigCount; e++)
    {
        if(form_factors[e] != NULL)
            free(form_factors[e]);
        if( e < AtomTypeCount && vacuo_form_factors[e] != NULL)
            free(vacuo_form_factors[e]);
    }
}


void unfitted_saxs(double* restrict intensities, Atom* atoms, size_t num_atoms,
                   double* restrict q_grid, size_t len_grid, coord* sampling_vectors, size_t num_vectors,
                   double solvent_density, double radius_adjustment)
{
    /*
     * Step 0: Pre-calculate any values that would benefit being taken out of the main nested loops, i.e. any components
     * of the equation that don't need to know q, v, and a simultaneously
     */
    
    // The form factors for each atom in a vacuum as depend only on q and atom.element
    double* vacuo_form_factors[AtomTypeCount] = { NULL };
    double* form_factors[AtomTypeCount + OrganicAtomConfigCount] = { NULL };
    size_t* offsets = malloc(sizeof(size_t) * num_atoms);   // maps indices of atom --> form_factors
    // Avoid repeated squaring of q when adjusting dummy atom form factors
    double* q_squared = malloc(sizeof(double) * len_grid);
    
    // These properties can be moved out of the q loop for calculating the corrected form factors
    double vol;
    double pivol23;
    double dummy_ff;
    
    
    clock_t begin, end;
    double time_spent;
    begin = clock();
    
    // do hydrogens first to add implicitly to combined organic atoms
    //      also populate the q squared lookup table while we're at it
    vacuo_form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    vol = sphere_volume(VdwRadii[Hydrogen] * radius_adjustment);
    pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
    dummy_ff = solvent_density * vol;
    for(size_t q = 0; q < len_grid; q++)
    {
        double q2 = q_grid[q] * q_grid[q];
        double vff = form_factor(Hydrogen, q_grid[q]);
        double ff = vff - dummy_ff * exp(q2 * pivol23);
        vacuo_form_factors[Hydrogen][q] = vff;
        form_factors[Hydrogen][q] = ff;
        q_squared[q] = q2;
    }
    
    for(size_t a = 0; a < num_atoms; a++)
    {
        offsets[a] = atoms[a].organic_type < UnknownBonds ? AtomTypeCount+atoms[a].organic_type : atoms[a].element;
        assert(offsets[a] >= AtomTypeCount);
        // If this is the first occurance of this element
        // ... then we can calculate the vacuo form factors across all q, while also doing the form factors for this atom
        if(vacuo_form_factors[atoms[a].element] == NULL)
        {
            
            vacuo_form_factors[atoms[a].element] = malloc(sizeof(double) * len_grid);
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            //vol = get_atom_volume(atoms+a,true);
            vol = get_atom_volume(atoms + a, true) * radius_adjustment * radius_adjustment * radius_adjustment;
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = form_factor(atoms[a].element, q_grid[q]);
                vacuo_form_factors[atoms[a].element][q] = vff;
                double ff = vff + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q] - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
        //  Else if we've never seen this combined organic atom before
        //  ... just do the corrected form factor calculation using previously stored vacuo factors
        else if(form_factors[offsets[a]] == NULL)
        {
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            vol = get_atom_volume(atoms + a, true) * radius_adjustment * radius_adjustment * radius_adjustment;
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = vacuo_form_factors[atoms[a].element][q] + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q];
                double ff = vff - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
    }
    // Step 1: Place solvent globs
    // Step 1: Place solvent globs
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Populating form factor table took %0.4f sec\n", time_spent);
    
    begin = clock();
    // Step 2: calculate scattering
    #pragma omp parallel for
    for(size_t q = 0; q < len_grid; q++)
    {
        double intensity = 0.0;
        for(size_t v = 0; v < num_vectors; v++)
        {
            double cos_bucket = 0.0;
            double sin_bucket = 0.0;
            coord q_vector = vec_multiply(sampling_vectors[v], q_grid[q]);
            
            for(size_t a = 0; a < num_atoms; a++)
            {
                Atom* atom = atoms+a;
                double corrected_f = form_factors[offsets[a]][q];
                double scat = dot_product(q_vector,atom->position);
                cos_bucket += corrected_f * cos(scat);
                sin_bucket += corrected_f * sin(scat);
            }
            intensity += cos_bucket * cos_bucket + sin_bucket * sin_bucket;
        }
        intensities[q] = intensity / num_vectors;
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Calculating scattering took %0.4f sec\n", time_spent);
    
    // Free allocated memory
    free(offsets);
    free(q_squared);
    for(size_t e = 0; e < AtomTypeCount + OrganicAtomConfigCount; e++)
    {
        if(form_factors[e] != NULL)
            free(form_factors[e]);
        if( e < AtomTypeCount && vacuo_form_factors[e] != NULL)
            free(vacuo_form_factors[e]);
    }
}

void debeye(double* restrict intensities, Atom* atoms, size_t num_atoms,
            double* restrict q_grid, size_t len_grid, double solvent_density)
{
    /*
    * Step 0: Pre-calculate any values that would benefit being taken out of the main nested loops, i.e. any components
    * of the equation that don't need to know q, v, and a simultaneously
    */
    // The form factors for each atom in a vacuum as depend only on q and atom.element
    double* vacuo_form_factors[AtomTypeCount] = { NULL };
    double* form_factors[AtomTypeCount + OrganicAtomConfigCount] = { NULL };
    size_t* offsets = malloc(sizeof(size_t) * num_atoms);   // maps indices of atom --> form_factors
    // Avoid repeated squaring of q when adjusting dummy atom form factors
    double* q_squared = malloc(sizeof(double) * len_grid);
    
    // These properties can be moved out of the q loop for calculating the corrected form factors
    double vol;
    double pivol23;
    double dummy_ff;
    
    // do hydrogens first to add implicitly to combined organic atoms
    //      also populate the q squared lookup table while we're at it
    vacuo_form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    form_factors[Hydrogen] = malloc(sizeof(double) * len_grid);
    vol = sphere_volume(VdwRadii[Hydrogen]);
    pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
    dummy_ff = solvent_density * vol;
    for(size_t q = 0; q < len_grid; q++)
    {
        double q2 = q_grid[q] * q_grid[q];
        double vff = form_factor(Hydrogen, q_grid[q]);
        double ff = vff - dummy_ff * exp(q2 * pivol23);
        vacuo_form_factors[Hydrogen][q] = vff;
        form_factors[Hydrogen][q] = ff;
        q_squared[q] = q2;
    }
    
    for(size_t a = 0; a < num_atoms; a++)
    {
        offsets[a] = atoms[a].organic_type < UnknownBonds ? AtomTypeCount+atoms[a].organic_type : atoms[a].element;
        // If this is the first occurance of this element
        // ... then we can calculate the vacuo form factors across all q, while also doing the form factors for this atom
        if(vacuo_form_factors[atoms[a].element] == NULL)
        {
            
            vacuo_form_factors[atoms[a].element] = malloc(sizeof(double) * len_grid);
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            vol = get_atom_volume(atoms+a,true);
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = form_factor(atoms[a].element, q_grid[q]);
                vacuo_form_factors[atoms[a].element][q] = vff;
                double ff = vff + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q] - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
            //  Else if we've never seen this combined organic atom before
            //  ... just do the corrected form factor calculation using previously stored vacuo factors
        else if(form_factors[offsets[a]] == NULL)
        {
            form_factors[offsets[a]] = malloc(sizeof(double) * len_grid);
            vol = get_atom_volume(atoms+a,true);
            pivol23 = -pow(vol, 2.0/3.0) / (coord_pi * 16.0);
            dummy_ff = solvent_density * vol;
            for(size_t q = 0; q < len_grid; q++)
            {
                double vff = vacuo_form_factors[atoms[a].element][q] + OrganicAtomHCount[atoms[a].organic_type] * vacuo_form_factors[Hydrogen][q];
                double ff = vff - dummy_ff * exp(q_squared[q] * pivol23);
                form_factors[offsets[a]][q] = ff;
            }
        }
    }
    
    #pragma omp parallel for
    for(size_t q = 0; q < len_grid; q++)
    {
        double intensity = 0.0;
        for(size_t a1 = 0; a1 < num_atoms; a1++)
        {
            double ff1 = form_factors[offsets[a1]][q];
            for(size_t a2 = 0; a2 < num_atoms; a2++)
            {
                if(a1 == a2)
                    intensity += ff1 * ff1;
                else
                {
                    double ff2 = form_factors[offsets[a2]][q];
                    if(q_grid[q] == 0.0)
                        intensity += ff1 * ff2;
                    else
                    {
                        double qd = distance(atoms[a1].position, atoms[a2].position) * q_grid[q];
                        intensity += ff1 * ff2 * sin(qd) / qd;
                    }
                }
            }
        }
        intensities[q] = intensity;
    }
}
