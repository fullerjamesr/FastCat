#pragma once
#ifndef FASTCAT_CPDB_H
#define FASTCAT_CPDB_H

#include <stdio.h>
#include <stdint.h>
#include "coord.h"
#include "atominfo.h"

static const char CPDB_DEFAULT_STR = '\0';

typedef enum ResType
{
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLU,
    GLN,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL,
    DA,
    DT,
    DG,
    DC,
    DU,
    RA,
    RT,
    RG,
    RC,
    RU,
    HOH,
    UnknownRes,
    
    ResTypeCount
} ResType;
extern const char* ResidueTypeStrings[];

typedef enum OrganicAtomConfig
{
    C3H0,
    C3H1,
    C4H1,
    C4H2,
    C4H3,
    N2H0,
    N3H0,
    N3HImidazole,
    N3H1,
    N3HGuanidinium,
    N3H2,
    N4H3,
    O1H0,
 //   O1H2, //WATER
    O2H0,
    O2H1,
    S2H0,
    S2H1,
    P4H0,
    
    UnknownBonds,
    OrganicAtomConfigCount
} OrganicAtomConfig;
extern const double OrganicAtomVdw[];
extern const double OrganicAtomVolumes[];
extern const unsigned int OrganicAtomHCount[];
extern const AtomType OrganicAtomBaseElement[];

typedef struct Atom
{
    uint32_t serial;            // atom serial number from PDB file
    char name[5];               // name (e.g. N, CA, CB, CG, etc... generally the element plus a specifier
    char res_name[4];           // the name of the residue that this atom belongs to (e.g. MET)
    char chain;                 // the one letter chainID this atom belongs to
    int res_num;                // the residue number this atom belongs to
    coord position;             // x,y,z coordinates
    float occupancy;            // from pdb file
    float temp;                 // B factor from pdb file
    char elementstr[3];         // character representation of element
    int8_t charge;              // charge on this atom
    
    AtomType element;           // better representation of element
    ResType res_type;           // the ResidueType of the parent of this Atom
    OrganicAtomConfig organic_type;
} Atom;



size_t read_atoms_from_pdbfile(Atom** buffer, size_t buffer_size, FILE* source, bool incl_hydrogens);

typedef struct NeighborList
{
    coord origin;               // The coordinates of the corner of the overall box that encloses all the atoms
    double edge_length;         // The edge length of the individual cells
    size_t dimensions[3];       // The number of cells [X, Y, Z] starting from `origin`
    size_t cell_lists_len;      // The total number of cells (X * Y * Z)
    size_t** cell_lists;        // For each cell, a list of the indices of atoms that belong to it
    size_t* cells_len;          // The length of each list in `cell lists`
    size_t* cell_assignments;   // `cell_assignments[a]` is the index in `cell lists` in which `atoms[a]` resides
} NeighborList;
size_t allocate_make_neighbor_list(NeighborList *nl, const Atom *atoms, const size_t num_atoms, const double edge_length);
void free_neighborlist(NeighborList* nl);

OrganicAtomConfig determine_bonding_config(const Atom* atom);
void assign_organic_atom_types(Atom* atoms, size_t num_atoms);
double get_vdw_radius(Atom* atom, bool implicit_hydrogen);
double get_atom_volume(Atom* atom, bool implicit_hydrogen);


double shrake_rupley_sasa(double* restrict sasa_per_atom, const double probe_radius,
                          const Atom* restrict atoms, const size_t num_atoms, bool implicit_hydrogen,
                          NeighborList* neighborList, coord* restrict sampling_vectors,
                          const size_t num_vectors);

#endif //FASTCAT_CPDB_H