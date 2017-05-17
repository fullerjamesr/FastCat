#pragma once
#ifndef FASTCAT_SAXS_H
#define FASTCAT_SAXS_H

#include "cpdb.h"

typedef struct SaxsProfile
{
    double* q;          // Q values
    double* restrict i; // I values
    double* restrict e; // Errors, can be null for theoretical profiles
    
    size_t length;      // The length of the q, i, e vectors
} SaxsProfile;

typedef struct ProfileScaleResult
{
    double constant;
    double offset;
} ProfileScaleResult;

typedef struct FormFactorTable
{
    double* restrict vacuo_form_factors;
    double* restrict dummy_form_factors;
    size_t dummy_options_count;
    double* restrict solvent_form_factors;
} FormFactorTable;

ProfileScaleResult find_chisq_scale_factor(const SaxsProfile *const restrict reference,
                                           const SaxsProfile *const restrict to_move);
ProfileScaleResult find_chisq_scale_factor_with_offset(const SaxsProfile *const restrict reference,
                                                       const SaxsProfile *const restrict to_move);

double chi_square(const SaxsProfile* const restrict experimental, const SaxsProfile* const restrict model);

double chi_square_fit(const SaxsProfile* const restrict experimental, const SaxsProfile* const restrict model);

double fast_form_factor(unsigned int z, double q);
double form_factor(unsigned int z, double q);

void allocate_form_factor_table(FormFactorTable* const to_init, const size_t q_length, const size_t radius_samples);

void populate_form_factor_table(FormFactorTable* const destination, const Atom* const atoms, const size_t atoms_length,
                                const double* const restrict q_values, const size_t q_values_length,
                                const double min_radius_adjustment, const double max_radius_adjustment, const size_t num_radius_samples,
                                const double dummy_solvent_density);


void unfitted_saxs(SaxsProfile* restrict saxs, const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                   const coord* restrict sampling_vectors, const size_t num_vectors,
                   const FormFactorTable* const formFactorTable, const size_t dummy_offset, const double solvation_factor);

void fitted_saxs(double* const restrict intensities, const double* const restrict q_values, const size_t profile_length,
                 const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                 const coord* restrict sampling_vectors, const size_t num_vectors,
                 const FormFactorTable* const formFactorTable, const double min_hydration, const double max_hydration, const size_t num_hydration_samples);

void unfitted_debeye(SaxsProfile* restrict saxs, const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                     const FormFactorTable* const formFactorTable, const size_t dummy_offset, const double solvation_factor);

#endif //FASTCAT_SAXS_H