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
    double q_spacing;   // The spacing between Q values
} SaxsProfile;

typedef struct FormFactorTable
{
    double* restrict vacuo_form_factors;
    double* restrict dummy_form_factors;
    size_t dummy_options_count;
    double* restrict solvent_form_factors;
} FormFactorTable;


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

void fitted_saxs(double* restrict intensities, const double* restrict q_grid, const size_t len_grid,
                   const Atom* restrict atoms, const double* restrict per_atom_sasa, const size_t num_atoms,
                   const coord* restrict sampling_vectors, const size_t num_vectors,
                   const double solvent_density, const double radius_adjustment);

void unfitted_debeye(SaxsProfile* restrict saxs, const Atom* const restrict atoms, const double* const restrict per_atom_sasa, const size_t num_atoms,
                     const FormFactorTable* const formFactorTable, const size_t dummy_offset, const double solvation_factor);

#endif //FASTCAT_SAXS_H