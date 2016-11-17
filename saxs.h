#pragma once
#ifndef FASTCAT_SAXS_H
#define FASTCAT_SAXS_H

#include "cpdb.h"

typedef struct SaxsProfile
{
    double* q;
    double* i;
    double* e;
    
    size_t length;
    double q_grid;
} SaxsProfile;

double fast_form_factor(unsigned int z, double q);
double form_factor(unsigned int z, double q);

void unfitted_saxs(double* restrict intensities, Atom* atoms, size_t num_atoms,
                   double* restrict q_grid, size_t len_grid, coord* sampling_vectors, size_t num_vectors,
                   double solvent_density, double radius_adjustment);

void unfitted_saxs2(double* restrict intensities, Atom* atoms, size_t num_atoms,
                   double* restrict q_grid, size_t len_grid, coord* sampling_vectors, size_t num_vectors,
                   double solvent_density, double radius_adjustment);

void debeye(double* restrict intensities, Atom* atoms, size_t num_atoms,
            double* restrict q_grid, size_t len_grid, double solvent_density);

#endif //FASTCAT_SAXS_H