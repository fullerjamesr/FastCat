#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "saxs.h"


static const double bohr = 0.52917721067;

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
//TODO: Having to increment z because Hydrogen = 0 will almost certainly be a source of confusion later
//TODO: Have a think about how to speed this calculation up further
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


// TODO: is size_t, which is generally 8 bytes, unnecessarily large for enums, which are generally 4 bytes?
static const size_t MAX_OFFSET = (size_t) (AtomTypeCount + OrganicAtomConfigCount);
static inline size_t atom2offset(const Atom *atom)
{
    return atom->organic_type < UnknownBonds ? AtomTypeCount + atom->organic_type : atom->element;
}

void allocate_form_factor_table(FormFactorTable* const to_init, const size_t q_length, const size_t radius_samples)
{
    assert(to_init != NULL);
    
    to_init->vacuo_form_factors = malloc(sizeof(double) * q_length * MAX_OFFSET);
    to_init->dummy_form_factors = malloc(sizeof(double) * q_length * MAX_OFFSET * radius_samples);
    to_init->solvent_form_factors = malloc(sizeof(double) * q_length * MAX_OFFSET);
    
    assert(to_init->vacuo_form_factors != NULL);
    assert(to_init->dummy_form_factors != NULL);
    assert(to_init->solvent_form_factors != NULL);
    
    to_init->dummy_options_count = radius_samples;
}
static inline void enumerate_dummy_volume_adjustments(double* const restrict destination, const double q,
                                                      const double default_volume, const double min_radius_adjustment,
                                                      const double max_radius_adjustment, const size_t num_radius_samples,
                                                      const double dummy_solvent_density)
{
    double delta = (max_radius_adjustment - min_radius_adjustment) / num_radius_samples;
    
    double constant = q * q / (coord_pi * 16.0);
    
    for(size_t i = 0; i < num_radius_samples; i++)
    {
        double adjustment = min_radius_adjustment + delta * i;
        double current_volume = default_volume * adjustment * adjustment * adjustment;
        destination[i] = dummy_solvent_density * current_volume * exp(-pow(current_volume, 2.0/3.0) * constant);
    }
}
/*
 *  Assumes
 *      destination->vacuo_form_factors is at least `q_values_length` * (AtomTypeCount + OrganicAtomTypeCount) long
 *      destination->dummy_form_factors is at least `q_values_length` * (AtomTypeCount + OrganicAtomTypeCount) * `num_radius_samples` long
 *      destination->solvent_form_factors is at least `q_values_length` long
 */
void populate_form_factor_table(FormFactorTable* const destination, const Atom* const atoms, const size_t atoms_length,
                                const double* const restrict q_values, const size_t q_values_length,
                                const double min_radius_adjustment, const double max_radius_adjustment, const size_t num_radius_samples,
                                const double dummy_solvent_density)
{
    // Holds a boolean list for identifying the set of elements and combined atoms present in `atoms`
    bool* restrict elements_to_build = calloc(AtomTypeCount, sizeof(bool));
    bool* restrict heavy_atoms_to_build = calloc(OrganicAtomConfigCount, sizeof(bool));
    // This boolean array is used in conjunction with needs_to_be_built, the difference being that base elements are
    // NOT set to true in this list, unlike `needs_to_be_built`
    bool* restrict exists_in_molecule = calloc(AtomTypeCount, sizeof(bool));
    
    for(size_t a = 0; a < atoms_length; a++)
    {
        elements_to_build[atoms[a].element] = true;
        if(atoms[a].organic_type < UnknownBonds)
            heavy_atoms_to_build[atoms[a].organic_type] = true;
        else
            exists_in_molecule[atoms[a].element] = true;
    }
    // Hydrogen is unlikely to be in a PDB but is needed for basically all heavy atoms
    elements_to_build[Hydrogen] = true;
    // Needed for water (for solvent)
    elements_to_build[Oxygen] = true;
    
    // Fill `destination->vacuo_form_factors` and `destination->dummy_form_factors`, but only for those elements/heavy
    // atoms present in `atoms`. The solvent (H2O) can also be done based on Oxygen and Hydrogen.
    for(size_t q = 0; q < q_values_length; q++)
    {
        size_t heavy_row = q * MAX_OFFSET;
        double q_value = q_values[q];
        
        for(size_t e = 0; e < AtomTypeCount; e++)
        {
            if(elements_to_build[e])
            {
                destination->vacuo_form_factors[heavy_row + e] = form_factor(e, q_value);
                if(exists_in_molecule[e])
                {
                    enumerate_dummy_volume_adjustments(destination->dummy_form_factors +
                                                       (q * MAX_OFFSET * num_radius_samples + e * num_radius_samples),
                                                       q_value, VdwSphereVolumes[e], min_radius_adjustment,
                                                       max_radius_adjustment, num_radius_samples, dummy_solvent_density);
                }
            }
        }
    
        for(size_t e = 0; e < OrganicAtomConfigCount; e++)
        {
            if(heavy_atoms_to_build[e])
            {
                size_t offset = e + AtomTypeCount;
                destination->vacuo_form_factors[heavy_row + offset] =
                        destination->vacuo_form_factors[heavy_row + OrganicAtomBaseElement[e]] +
                        destination->vacuo_form_factors[heavy_row + Hydrogen] * OrganicAtomHCount[e];
                enumerate_dummy_volume_adjustments(destination->dummy_form_factors +
                                                       (heavy_row * num_radius_samples + offset * num_radius_samples),
                                                       q_value, OrganicAtomVolumes[e], min_radius_adjustment,
                                                       max_radius_adjustment, num_radius_samples, dummy_solvent_density);
            }
        }
        
        // The solvent goes below. The assumption is it's H2O, and the constants below are pre-calculated assuming a
        // radius of 1.525 angstroms.
        destination->solvent_form_factors[q] = destination->vacuo_form_factors[heavy_row + Oxygen]
                                               + destination->vacuo_form_factors[heavy_row + Hydrogen] * 2
                                               - dummy_solvent_density * 14.85587171 * exp(q_value * q_value * 0.1202252175);
    }

    free(elements_to_build);
    free(heavy_atoms_to_build);
    free(exists_in_molecule);
}


//TODO: Consider moving the spooling out of positions and offsets to a separate function. This could also contain the form factor table.
void unfitted_saxs(SaxsProfile* restrict saxs, const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                   const coord* restrict sampling_vectors, const size_t num_vectors,
                   const FormFactorTable* const formFactorTable, const size_t dummy_offset, const double solvation_factor)
{
    const size_t saxs_profile_length = saxs->length;
    const size_t dummy_options_count = formFactorTable->dummy_options_count;
 
    // spool out offsets and positions for better memory layout
    coord* restrict positions = malloc(sizeof(coord) * num_atoms);
    size_t* restrict offsets = malloc(sizeof(size_t) * num_atoms);
    for(size_t a = 0; a < num_atoms; a++)
    {
        positions[a] = atoms[a].position;
        offsets[a] = atom2offset(atoms+a);
    }
    
    #pragma omp parallel for
    for(size_t q = 0; q < saxs_profile_length; q++)
    {
        double intensity = 0.0;
        double q_value = saxs->q[q];
        size_t ff_table_row = q * MAX_OFFSET;
        for(size_t v = 0; v < num_vectors; v++)
        {
            double cos_bucket = 0.0;
            double sin_bucket = 0.0;
            coord q_vector = vec_multiply(sampling_vectors[v], q_value);
            
            for(size_t a = 0; a < num_atoms; a++)
            {
                double corrected_f = formFactorTable->vacuo_form_factors[ff_table_row + offsets[a]]
                                     - formFactorTable->dummy_form_factors[ff_table_row * dummy_options_count + offsets[a] * dummy_options_count + dummy_offset]
                                     + per_atom_sasa[a] * solvation_factor * formFactorTable->solvent_form_factors[q];
                double scat = dot_product(q_vector, positions[a]);
                //TODO: Using the trig identity to get sin is faster on both gcc and icc, but we should find out if this
                //      is still the case after accounting for the fact that we need to keep track of quadrant
                double cos_scat = cos(scat);
                double sin_scat = sin(scat);
                cos_bucket += corrected_f * cos_scat;
                sin_bucket += corrected_f * sin_scat;
            }
            intensity += cos_bucket * cos_bucket + sin_bucket * sin_bucket;
        }
        saxs->i[q] = intensity / num_vectors;
    }
    
    free(offsets);
    free(positions);
}

void fitted_saxs(double* const restrict intensities, const double* const restrict q_values, const size_t profile_length,
                 const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                 const coord* restrict sampling_vectors, const size_t num_vectors,
                 const FormFactorTable* const formFactorTable, const double min_hydration, const double max_hydration, const size_t num_hydration_samples)
{
    // Sanity input checks for debugging
    assert(intensities != NULL);
    assert(q_values != NULL);
    assert(atoms != NULL);
    assert(per_atom_sasa != NULL);
    assert(sampling_vectors != NULL);
    assert(formFactorTable != NULL);
    assert(num_hydration_samples > 0);
    assert(num_atoms > 0);
    
    const size_t dummy_options_count = formFactorTable->dummy_options_count;
    const double delta_hydration = (max_hydration - min_hydration) /num_hydration_samples;
    
    // spool out offsets and positions for better memory layout
    coord* const restrict positions = malloc(sizeof(coord) * num_atoms);
    size_t* const restrict offsets = malloc(sizeof(size_t) * num_atoms);
    for(size_t a = 0; a < num_atoms; a++)
    {
        positions[a] = atoms[a].position;
        offsets[a] = atom2offset(atoms+a);
    }

    double* restrict intensity = calloc(dummy_options_count * num_hydration_samples, sizeof(double));
    double* restrict dummy_cos_bucket = calloc(dummy_options_count, sizeof(double));
    double* restrict dummy_sin_bucket = calloc(dummy_options_count, sizeof(double));
    assert(intensity[0] == 0.0);
    assert(dummy_cos_bucket[0] == 0.0);
    assert(dummy_sin_bucket[0] == 0.0);
    
    #pragma omp parallel for
    for(size_t q = 0; q < profile_length; q++)
    {
        size_t ff_table_row = q * MAX_OFFSET;
        double* restrict vacuo_form_factors = formFactorTable->vacuo_form_factors + ff_table_row;
        for(size_t v = 0; v < num_vectors; v++)
        {
            coord q_vector = vec_multiply(sampling_vectors[v], q_values[q]);
            
            double vacuum_cos_bucket = 0.0;
            double vacuum_sin_bucket = 0.0;
            double hydration_cos_bucket = 0.0;
            double hydration_sin_bucket = 0.0;
            
            for(size_t a = 0; a < num_atoms; a++)
            {
                double scat = dot_product(q_vector, positions[a]);
                double cos_scat = cos(scat);
                double sin_scat = sin(scat);
    
                double vacuum_ff = vacuo_form_factors[ offsets[a] ];
                vacuum_cos_bucket += vacuum_ff * cos_scat;
                vacuum_sin_bucket += vacuum_ff * sin_scat;
    
                double hydration_ff = per_atom_sasa[a] * formFactorTable->solvent_form_factors[q];
                hydration_cos_bucket += hydration_ff * cos_scat;
                hydration_sin_bucket += hydration_ff * sin_scat;
                
                double* restrict dummy_ffs = formFactorTable->dummy_form_factors + ff_table_row * dummy_options_count + offsets[a] * dummy_options_count;
                for(size_t d = 0; d < dummy_options_count; d++)
                {
                    double dummy_ff = dummy_ffs[d];
                    dummy_cos_bucket[d] += dummy_ff * cos_scat;
                    dummy_sin_bucket[d] += dummy_ff * sin_scat;
                }
            }
            
            for(size_t d = 0; d < dummy_options_count; d++)
            {
                for(size_t h = 0; h < num_hydration_samples; h++)
                {
                    double cosine_term = vacuum_cos_bucket - dummy_cos_bucket[d] +
                                         hydration_cos_bucket * (min_hydration + delta_hydration * h);
                    double sin_term = vacuum_sin_bucket - dummy_sin_bucket[d] +
                                      hydration_sin_bucket * (min_hydration + delta_hydration * h);
                    intensity[d * num_hydration_samples + h] += cosine_term * cosine_term + sin_term * sin_term;
                }
            }
    
            memset(dummy_cos_bucket, 0, sizeof(double) * dummy_options_count);
            memset(dummy_sin_bucket, 0, sizeof(double) * dummy_options_count);
            assert(dummy_cos_bucket[0] == 0.0);
            assert(dummy_sin_bucket[0] == 0.0);
        }
        
        // Spool the results of the calculation back into `intensities`, the caller-supplied destination array
        // The memory access pattern here is bad in order to keep the memory access pattern of the previous loops good,
        // because dummy_options_count * num_hydration_samples should be <<<< num_vectors * num_atoms
        for(size_t d = 0; d < dummy_options_count; d++)
        {
            size_t inner_dim_offset = d * num_hydration_samples;
            size_t outer_dim_offset = inner_dim_offset * profile_length + q;
            for(size_t h = 0; h < num_hydration_samples; h++)
            {
                intensities[outer_dim_offset + h * profile_length] = intensity[inner_dim_offset + h] / num_vectors;
            }
        }
    
        memset(intensity, 0, sizeof(double) * dummy_options_count * num_hydration_samples);
        assert(intensity[0] == 0.0);
    }
    
    free(offsets);
    free(positions);
    free(intensity);
    free(dummy_cos_bucket);
    free(dummy_sin_bucket);
}

void unfitted_debeye(SaxsProfile* restrict saxs, const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                     const FormFactorTable* const formFactorTable, const size_t dummy_offset, const double solvation_factor)
{
    const size_t saxs_profile_length = saxs->length;
    const size_t dummy_options_count = formFactorTable->dummy_options_count;
    
    // spool out offsets and positions for better memory layout
    coord* restrict positions = malloc(sizeof(coord) * num_atoms);
    size_t* restrict offsets = malloc(sizeof(size_t) * num_atoms);
    for(size_t a = 0; a < num_atoms; a++)
    {
        positions[a] = atoms[a].position;
        offsets[a] = atom2offset(atoms+a);
    }
    
    #pragma omp parallel for
    for(size_t q = 0; q < saxs_profile_length; q++)
    {
        double intensity = 0.0;
        double q_value = saxs->q[q];
        size_t ff_table_row = q * MAX_OFFSET;
        for(size_t a1 = 0; a1 < num_atoms; a1++)
        {
            double ff1 = formFactorTable->vacuo_form_factors[ff_table_row + offsets[a1]]
                         - formFactorTable->dummy_form_factors[ff_table_row * dummy_options_count + offsets[a1] * dummy_options_count + dummy_offset]
                         + per_atom_sasa[a1] * solvation_factor * formFactorTable->solvent_form_factors[q];
            for(size_t a2 = 0; a2 < num_atoms; a2++)
            {
                if(a1 == a2)
                    intensity += ff1 * ff1;
                else
                {
                    double ff2 = formFactorTable->vacuo_form_factors[ff_table_row + offsets[a1]]
                                 - formFactorTable->dummy_form_factors[ff_table_row * dummy_options_count + offsets[a1] * dummy_options_count + dummy_offset]
                                 + per_atom_sasa[a1] * solvation_factor * formFactorTable->solvent_form_factors[q];
                    // Avoid problems with evaluation of the sinc function if q is approaching 0
                    if(q_value < 0.0001)
                        intensity += ff1 * ff2;
                    else
                    {
                        double qd = distance(atoms[a1].position, atoms[a2].position) * q_value;
                        intensity += ff1 * ff2 * sin(qd) / qd;
                    }
                }
            }
        }
        saxs->i[q] = intensity;
    }
    
    free(offsets);
    free(positions);
}
