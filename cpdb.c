#include <assert.h>
#include "fileio.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "cpdb.h"

const char* ResidueTypeStrings[] =
{
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "DA",
    "DT",
    "DG",
    "DC",
    "DU",
    "A",
    "T",
    "G",
    "C",
    "U",
    "HOH"
};

const double OrganicAtomVdw[] =
{
    1.61,   // C3H0
    1.76,   // C3H1
    1.88,   // C4H1
    1.88,   // C4H2
    1.88,   // C4H3
    1.61,   // N2H0
    1.64,   // N3H0
    1.64,   // N3HImidazole
    1.64,   // N3H1
    1.64,   // N3HGuanidinium
    1.64,   // N3H2
    1.64,   // N4H3
    1.42,   // O1H0
    1.51,   // O2H0
    1.46,   // O2H1
    1.77,   // S2H0
    1.77,   // S2H1
    2.04    // P4H0
};
const AtomType OrganicAtomBaseElement[] =
{
    Carbon,     // C3H0
    Carbon,     // C3H1
    Carbon,     // C4H1
    Carbon,     // C4H2
    Carbon,     // C4H3
    Nitrogen,   // N2H0
    Nitrogen,   // N3H0
    Nitrogen,   // N3HImidazole
    Nitrogen,   // N3H1
    Nitrogen,   // N3HGuanidinium
    Nitrogen,   // N3H2
    Nitrogen,   // N4H3
    Oxygen,     // O1H0
    Oxygen,     // O2H0
    Oxygen,     // O2H1
    Sulfur,     // S2H0
    Sulfur,     // S2H1
    Phosphorus  // P4H0
};
const unsigned int OrganicAtomHCount[] =
{
    0,  // C3H0
    1,  // C3H1
    1,  // C4H1
    2,  // C4H2
    3,  // C4H3
    0,  // N2H0
    0,  // N3H0
    1,  // N3HImidazole
    1,  // N3H1
    2,  // N3HGuanidinium
    2,  // N3H2
    3,  // N4H3
    0,  // O1H0
    0,  // O2H0
    1,  // O2H1
    0,  // S2H0
    1,  // S2H1
    0   // P4H0
};
const double OrganicAtomVolumes[] =
{
    17.48099857,    // C3H0
    22.83634591,    // C3H1
    27.83313699,    // C4H1
    27.83313699,    // C4H2
    27.83313699,    // C4H3
    17.48099857,    // N2H0
    18.47651902,    // N3H0
    18.47651902,    // N3HImid
    18.47651902,    // N3H1
    18.47651902,    // N3HGuan
    18.47651902,    // N3H2
    18.47651902,    // N4H3
    11.99371273,    // O1H0
    14.42179942,    // O2H0
    13.03608479,    // O2H1
    23.22781767,    // S2H0
    23.22781767,    // S2H1
    35.56142141     // P4H0
};


/*
 * "Zero-out" an Atom struct, returning all member values to 0 / Unknown / empty values
 * This can also act as a constructor, if called with NULL
 */
Atom* initialize_atom(Atom* to_initialize)
{
    if(to_initialize == NULL)
    {
        to_initialize = malloc(sizeof(Atom));
        assert(to_initialize != NULL);
    }
    
    to_initialize->serial= 0;
    strcpy(to_initialize->name,&CPDB_DEFAULT_STR);
    strcpy(to_initialize->res_name,&CPDB_DEFAULT_STR);
    to_initialize->res_num = 0;
    to_initialize->chain= CPDB_DEFAULT_STR;
    to_initialize->position = UNITIALIZED_COORD;
    to_initialize->occupancy = 0.0f;
    to_initialize->temp = 0.0f;
    strcpy(to_initialize->elementstr, &CPDB_DEFAULT_STR);
    to_initialize->charge = 0;
    
    to_initialize->element = UnknownAtom;
    to_initialize->res_type = UnknownRes;
    to_initialize->organic_type = UnknownBonds;
    
    return to_initialize;
}

/*
 * Parse a line from a PDB file into an Atom struct.
 * Format is expected to conform to PDB spec http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
 * Does not determine information not immediately evident from the PDB text (e.g., covalent bonding configuration)
 *
 * Parameters:
 *  char* line: A properly null-terminated string that is one line conforming to the PDB format.
 *              Generally the line should start with ATOM or HETATM but this assumption is not checked
 *  Atom* atom: the struct to fill with the information parsed from `line`. It's contents will be re-initialized.
 *
 * Returns:
 *      true if there was at least enough data in the PDB line to determine an element and location from `line`
 *      false otherwise -- although there still may be some garbage values in `atom`
 */
bool parse_pdb_line_into_atom(const char* line, Atom * atom)
{
    // SANITY CHECKING
    assert(line != NULL);
    assert(atom != NULL);
    assert(strstr(line,"ATOM") == line || strstr(line,"HETATM") == line);
    // a proper pdb file is useless unless it at least has the x,y,z coords, which end at column 54
    size_t len_line=strlen(line);
    if(len_line < 54)
        return false;
    
    atom = initialize_atom(atom);
    
    // COLS 7-11 are the integer serial number
    atom->serial = (uint32_t) atoi(line + 6);
    // COLS 13-16 are the atom name
    strncpy(atom->name, line + 12, 4);
    atom->name[4] = '\0';
    // COLS 18-20 are the 3 letter residue codename
    strncpy(atom->res_name, line + 17, 3);
    atom->res_name[3]='\0';
    // COL 22 is the chain ID
    atom->chain=*(line + 21);
    // COLS 23-26 are the res number
    // need to copy this out because it abuts another field and I don't trust atoi
    char res_num[5];
    strncpy(res_num,line+22,4);
    res_num[4]='\0';
    atom->res_num=atoi(res_num);
    // COLS 31-38 are x coord, 39-46 y, 47-54
    // need to copy these out and null terminate before feeding to parsing functions because the spec
    // doesn't guarantee delimiting characters or whitespace
    char x[9],y[9],z[9];
    strncpy(x,line+30,8);
    strncpy(y,line+38,8);
    strncpy(z,line+46,8);
    x[8]='\0';y[8]='\0';z[8]='\0';
    char *xoff,*yoff,*zoff;
    atom->position.x = strtod(x, &xoff);
    if(xoff == x)
        return false;
    atom->position.y= strtod(y, &yoff);
    if(yoff == y)
        return false;
    atom->position.z= strtod(z, &zoff);
    if(zoff == z)
        return false;
    if(len_line > 59)
    {
        // COLS 55-60 and 61-66 are occupancy and temp factor and are similarly not guaranteed to be delimited
        char o[7];
        strncpy(o, line + 54, 6);
        o[6] = '\0';
        atom->occupancy = (float) atof(o);
        if (len_line > 65)
        {
            char t[7];
            strncpy(t, line + 60, 6);
            t[6] = '\0';
            atom->temp = (float) atof(t);
            if (len_line > 77)
            {
                // COLS 77-78 are the element symbol only (e.g. N, not NG), right justified
                strncpy(atom->elementstr, line + 76, 2);
                atom->elementstr[2]='\0';
                if (len_line > 79)
                {
                    // COLS 79-80 are the atoms charge (eg 1+, 2-)
                    char charge[2];
                    charge[0] = *(line + 78);
                    charge[1] = '\0';
                    atom->charge = (int8_t) atoi(charge);
                    if(*(line + 79) == '-')
                        atom->charge *= -1;
                }
            }
        }
    }
    
    // Determine the AtomType of this Atom
    // the most diagnostic would be current_atom->elementstr, if the line included it
    if(atom->elementstr[0] != '\0')
    {
        int j;
        for(j=0;j<(AtomTypeCount-1);j++)
        {
            // " if the AtomTypeString is in there and properly justified"
            // ... without the 2-len part below, cannot distinguish e.g. "C" from "CA"
            // convert to uppercase because disagreements appear to exist about whether it's OK to have lowercase here
            char element_toupper[3] = { (char)toupper(atom->elementstr[0]), (char)toupper(atom->elementstr[1]), '\0' };
            if(strstr(element_toupper, AtomicSymbolStrings[j]) == element_toupper + (2 - strlen(AtomicSymbolStrings[j])))
            {
                atom->element = (AtomType) j;
                break;
            }
        }
    }
    // ... but if needed it can be inferred from name
    else
    {
        // justification here is a nightmare of rules:
        // (1) if the first character is blank, then the element has a single letter symbol and is found at the second position
        // (2) if the last two characters are blank, then the element has a two letter symbol given by the first two positions
        // (3) finally, if neither of the above are true, the elements has a single letter symbol which is found at the first position
        int leading_space = isspace(atom->name[0]) ? 1 : 0;
        if(leading_space || (isdigit(atom->name[2]) && isdigit(atom->name[3])))
        {
            int i;
            for(i = 0; i < (AtomTypeCount - 1); i++)
            {
                // break on non-single-letter
                if(strlen(AtomicSymbolStrings[i]) > 1)
                    continue;
                if(atom->name[leading_space] == AtomicSymbolStrings[i][0])
                {
                    atom->element = (AtomType) i;
                    break;
                }
            }
    
            if(i >= (AtomTypeCount - 1))
                // there was no element specified and the atom name did not provide anything useful
                return false;
        }
        else
        {
            int i;
            for(i = AtomTypeCount - 2; i >= 0; i--)
            {
                if(strstr(atom->name, AtomicSymbolStrings[i]) == atom->name)
                {
                    atom->element=(AtomType) i;
                    break;
                }
            }
            
            if(i < 0)
                // there was no element specified and the atom name did not provide anything useful
                return false;
        }
    }
    
    // Determine the ResType of this atom
    for(int i = 0; i < ResTypeCount - 1; i++)
    {
        if (strstr(atom->res_name, ResidueTypeStrings[i]) != NULL)
        {
            atom->res_type = (ResType) i;
            break;
        }
    }
    
    return true;
}


/*
* Parse a list of Atoms in from a text PDB file.
*
* Parameters:
*    Atom** buffer: A pointer to a preallocated list of Atoms to read into
*    size_t buffer_size: The size of the preallocated list (in Atoms, NOT BYTES)
*    FILE* source: A FILE pointer open for reading to a PDB-formatted text file
*    bool incl_hydrogens: Whether or not to include hydrogens
*
* Returns:
*    The number of Atoms successfully read from `source` and saved in `buffer`
*/
size_t read_atoms_from_pdbfile(Atom** buffer, size_t buffer_size, FILE* source, bool incl_hydrogens)
{
    assert(buffer);
    assert(source);
    Atom* atoms = *buffer;
    assert(atoms);
    
    // The PDB specification says that any line in a PDB file is a maximum of 80 chars long, so this should be plenty
    char line[128] = { '\0' };
    size_t atoms_read = 0;
    while(!ferror(source) && !feof(source))
    {
        if(fgets(line, 128, source) && (strstr(line,"ATOM") == line || strstr(line,"HETATM") == line))
        {
            Atom* current_atom = atoms + atoms_read;
            if(parse_pdb_line_into_atom(line,current_atom) && (incl_hydrogens || current_atom->element != Hydrogen))
            {
                atoms_read++;
                if(atoms_read >= buffer_size)
                {
                    size_t newsize = buffer_size * 2;
                    Atom* newatoms = realloc(atoms,sizeof(Atom) * newsize);
                    if(newatoms == NULL)
                    {
                        // TODO: error message?
                        break;
                    }
                    atoms = newatoms;
                    buffer_size = newsize;
                }
            }
            // `parse_pdb_line_into_atom` will write as much info as it can get, so on rejection `current_atom` is wiped
            else
                initialize_atom(current_atom);
        }
    }
    
    *buffer = atoms;
    return atoms_read;
}


/*
 * Determine neighboring atoms by breaking the volume occupied by the structure into cubic cells
 *
 * Parameters:
 *  Atom* atoms: The list of atoms to consider
 *  size_t num_atoms: The number of atoms in `atoms` to consider
 *  double edge_length: The size of the cells to divide the atoms into
 *  size_t atom_cell_assignments: A preallocated array othat is at least as long as `num_atoms`. The index of the cell
 *      in which each atoms resides will be written here: atom_cell_assignments[i] = the cell in which atoms[i] exists
 *  size_t*** cells: A pointer to a 2D array to store the atoms the belong to each cell. cells[i] will be a list of
 *      indices in `atoms`. The value of this pointer will be allocated.
 *  size_t** cell_sizes: A pointer to an array. Will be allocated and the length of cells[i] will be recorded in
 *      cell_sizes[i].
 *
 * Returns:
 *  The total number of cells, which is the allocated size of `cells` and `cell_sizes`
 */
static inline size_t neighborlist_indexof(NeighborList *nl, size_t x, size_t y, size_t z)
{
    return x + (y * nl->dimensions[0]) + (z * nl->dimensions[0] * nl->dimensions[1]);
}
static inline void neighborlist_addressof(size_t* address, size_t idx, NeighborList* nl)
{
    // the z coordinate is idx / (dim[0] * dim[1])
    address[2] = idx / (nl->dimensions[0] * nl->dimensions[1]);
    // the y coordinate is the remainder of the above divided by dim[0]
    address[1] = (idx % (nl->dimensions[0] * nl->dimensions[1])) / nl->dimensions[0];
    // the x coordinate is the remainder of the above
    address[0] = (idx % (nl->dimensions[0] * nl->dimensions[1])) % nl->dimensions[0];
}
size_t allocate_make_neighbor_list(NeighborList *nl, const Atom *atoms, const size_t num_atoms,
                                   const double edge_length)
{
    // Sanity checks
    assert(atoms != NULL);
    assert(edge_length > 0.0);
    assert(nl != NULL);
    
    nl->edge_length = edge_length;
    
    // Step 1: Find the maximum dimensions along each axis to define a simple bounding cube
    double min_x = atoms[0].position.x;
    double max_x = min_x;
    double min_y = atoms[0].position.y;
    double max_y = min_y;
    double min_z = atoms[0].position.z;
    double max_z = min_z;
    for (size_t a = 1; a < num_atoms; a++)
    {
        if(atoms[a].position.x > max_x)
            max_x = atoms[a].position.x;
        else if(atoms[a].position.x < min_x)
            min_x = atoms[a].position.x;
        
        if(atoms[a].position.y > max_y)
            max_y = atoms[a].position.y;
        else if (atoms[a].position.y < min_y)
            min_y = atoms[a].position.y;
        
        if(atoms[a].position.z > max_z)
            max_z = atoms[a].position.z;
        else if(atoms[a].position.z < min_z)
            min_z = atoms[a].position.z;
    }
    
    nl->origin.x = min_x;
    nl->origin.y = min_y;
    nl->origin.z = min_z;
    
    // Step 2: Divide the bounding box into `edge_length`-sized cubic buckets
    size_t x_cells = (size_t) ceil((max_x - min_x) / edge_length);
    size_t y_cells = (size_t) ceil((max_y - min_y) / edge_length);
    size_t z_cells = (size_t) ceil((max_z - min_z) / edge_length);
    size_t num_cells = x_cells * y_cells * z_cells;
    
    nl->dimensions[0] = x_cells;
    nl->dimensions[1] = y_cells;
    nl->dimensions[2] = z_cells;
    nl->cell_lists_len = num_cells;
    
    // Step 3: Allocate buckets
    // `cell_lists` is a linear representation of the cube of buckets, i = x + (y * x_cells) + (z * x_cells * y_cells)
    // `cells_len` keeps track of how many atoms have been added to each `cell_bucket`
    nl->cell_lists = malloc(sizeof(size_t*) * num_cells);
    nl->cells_len = malloc(sizeof(size_t) * num_cells);
    nl->cell_assignments = malloc(sizeof(size_t) * num_atoms);
    for(size_t i = 0; i < num_cells; i++)
    {
        nl->cell_lists[i] = malloc(sizeof(size_t) * num_atoms);
        nl->cells_len[i] = 0;
    }
    
    //Step 4: Sort atoms into bucket by position
    for(size_t a = 0; a < num_atoms; a++)
    {
        size_t x_loc = (size_t) floor((atoms[a].position.x - min_x) / edge_length);
        size_t y_loc = (size_t) floor((atoms[a].position.y - min_y) / edge_length);
        size_t z_loc = (size_t) floor((atoms[a].position.z - min_z) / edge_length);
        size_t cell_address = neighborlist_indexof(nl, x_loc, y_loc, z_loc);
        nl->cell_assignments[a] = cell_address;
        nl->cell_lists[cell_address][nl->cells_len[cell_address]++] = a;
    }
    
    return num_cells;
}
size_t neighborlist_get_adjacent_cell_indices(size_t* dest, size_t idx, NeighborList* nl)
{
    assert(dest != NULL);
    assert(idx < nl->cell_lists_len);
    
    size_t address[3];
    neighborlist_addressof(address, idx, nl);
    
    int xmin = address[0] > 0 ? -1 : 0;
    int xmax = address[0] < (nl->dimensions[0] - 1) ? 1 : 0;
    int ymin = address[1] > 0 ? -1 : 0;
    int ymax = address[1] < (nl->dimensions[1] - 1) ? 1 : 0;
    int zmin = address[2] > 0 ? -1 : 0;
    int zmax = address[2] < (nl->dimensions[2] - 1) ? 1 : 0;
    size_t j = 0;
    for(int x = xmin; x < xmax; x++)
        for(int y = ymin; y < ymax; y++)
            for(int z = zmin; z < zmax; z++)
                if( ! (x == 0 && y == 0 && z == 0))
                    dest[j++] = neighborlist_indexof(nl, address[0] + x, address[1] + y, address[2] + z);
    
    return j;
}

void free_neighborlist(NeighborList* nl)
{
    for(int i = 0; i < nl->cell_lists_len; i++)
        free(nl->cell_lists[i]);
    free(nl->cell_lists);
    free(nl->cells_len);
    free(nl->cell_assignments);
}

OrganicAtomConfig determine_bonding_config(const Atom* atom)
{
    assert(atom != NULL);
    assert(atom->element == Carbon ||atom->element == Nitrogen ||
                   atom->element == Oxygen || atom->element == Sulfur || atom->element == Phosphorus);
    
    // TODO: I suppose it would be possible to brute force an answer by looking at nearby Atoms and applying some logic
    if(atom->res_type == UnknownRes)
        return UnknownBonds;
    
    // Peptide backbone
    if(strcmp(atom->name," CA ") == 0)
        return atom->res_type == GLY ? C4H2 : C4H1;
    else if(strcmp(atom->name," C  ") == 0)
        return C3H0;
    else if(strcmp(atom->name," O  ") == 0)
        return O1H0;
    else if(strcmp(atom->name," OXT") == 0)
        return O1H0;
    else if(strcmp(atom->name," OT1") == 0)
        return O1H0;
    else if(strcmp(atom->name," OT2") == 0)
        return O1H0;
    else if(strcmp(atom->name," N  ") == 0)
        return N3H1;
    //TODO: what is the most straightforward way to determine N-terminal-ness? Possibly related to the brute force above
    
    switch(atom->res_type)
    {
        case ALA:
            if(strcmp(atom->name," CB ") == 0)
                return C4H3;
            break;
        case ARG:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CD ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CZ ") == 0)
                return C3H0;
            else if (strcmp(atom->name, " NE ") == 0)
                return N3H1;
            else if(strcmp(atom->name, " NH1") == 0)
                return N3HGuanidinium;
            else if(strcmp(atom->name, " NH2") == 0)
                return N3HGuanidinium;
            break;
        case ASN:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " OD1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " ND2") == 0)
                return N3H2;
            break;
        case ASP:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " OD1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " OD2") == 0)
                return O1H0;
            break;
        case CYS:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " SG ") == 0)
                return S2H1;
            break;
        case GLU:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CD ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " OE1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " OE2") == 0)
                return O1H0;
            break;
        case GLN:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CD ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " OE1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " NE2") == 0)
                return N3H2;
            break;
        case HIS:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " ND1") == 0)
                return N3HImidazole;
            else if(strcmp(atom->name, " CD2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CE1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " NE2") == 0)
                return N3HImidazole;
            break;
        case ILE:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H1;
            else if(strcmp(atom->name, " CG1") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG2") == 0)
                return C4H3;
            else if(strcmp(atom->name, " CD1") == 0)
                return C4H3;
            break;
        case LEU:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H1;
            else if(strcmp(atom->name, " CD1") == 0)
                return C4H3;
            else if(strcmp(atom->name, " CD2") == 0)
                return C4H3;
        case LYS:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CD ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CE ") == 0)
                return C4H2;
            else if (strcmp(atom->name, " NZ ") == 0)
                return N4H3;
            break;
        case MET:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " SD ") == 0)
                return S2H0;
            else if(strcmp(atom->name, " CE ") == 0)
                return C4H3;
            break;
        case PHE:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " CD1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CD2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CE1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CE2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CZ ") == 0)
                return C3H1;
            break;
        case PRO:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CD ") == 0)
                return C4H2;
            break;
        case SER:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " OG ") == 0)
                return O2H1;
            break;
        case THR:
            if(strcmp(atom->name, " CB ") == 0)
                return C3H1;
            else if(strcmp(atom->name, " OG1") == 0)
                return O2H1;
            else if(strcmp(atom->name, " CG2") == 0)
                return C4H3;
            break;
        case TRP:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " CD1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CD2") == 0)
                return C3H0;
            else if(strcmp(atom->name, " NE1") == 0)
                return N3H1;
            else if(strcmp(atom->name, " CE2") == 0)
                return C3H0;
            else if(strcmp(atom->name, " CE3") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CZ2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CZ3") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CH2") == 0)
                return C3H1;
            break;
        case TYR:
            if(strcmp(atom->name, " CB ") == 0)
                return C4H2;
            else if(strcmp(atom->name, " CG ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " CD1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CD2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CE1") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CE2") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CZ ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " OH ") == 0)
                return O2H1;
            break;
        case VAL:
            if(strcmp(atom->name, " CB ") == 0)
                return C3H1;
            else if(strcmp(atom->name, " CG1") == 0)
                return C4H3;
            else if(strcmp(atom->name, " CG2") == 0)
                return C4H3;
            break;
        case RT:
        case DT:
            if(strcmp(atom->name, " C5 ") == 0)
                return C3H0;
        case RU:
        case DU:
            if(strcmp(atom->name, " N3 ") == 0)
                return N3H1;
            else if(strcmp(atom->name, " O4 ") == 0)
                return O1H0;
            else if(strcmp(atom->name, " C7 ") == 0)
                return C4H3;
        case DC:
        case RC:
            if(strcmp(atom->name, " N3 ") == 0)
                return N2H0;
            else if(strcmp(atom->name, " N4 ") == 0)
                return N3H2;
            else if(strcmp(atom->name, " C5 ") == 0)
                return C3H1;
            
            // Start common features of main purine ring
            if(strcmp(atom->name, " N1 ") == 0)
                return N3H0;
            else if(strcmp(atom->name, " C2 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " C4 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " C6 ") == 0)
                return C3H1;
            else if(strcmp(atom->name, " O2 ") == 0)
                return O1H0;
            // End common features of main purine ring
            // Start phosphate/sugar backbone
            if(strcmp(atom->name, " P  ") == 0)
                return P4H0;
            else if(strcmp(atom->name, " OP1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " OP2") == 0)
                return O1H0;
            else if(strcmp(atom->name, " O5'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " C5'") == 0)
                return C4H2;
            else if(strcmp(atom->name, " C4'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " C3'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " C2'") == 0)
                if(atom->res_type <= DU)
                    return C4H2;
                else
                    return C4H1;
            else if(strcmp(atom->name, " C1'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " O4'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " O3'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " O2'") == 0)
                return O2H1;
            break;
        // End phosphate/sugar backbone
        case DG:
        case RG:
            if(strcmp(atom->name, " O6 ") == 0)
                return O1H0;
            else if(strcmp(atom->name, " C2 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " N2 ") == 0)
                return N3H2;
            else if(strcmp(atom->name, " N1 ") == 0)
                return N3H1;
        case DA:
        case RA:
            if(strcmp(atom->name, " N6 ") == 0)
                return N3H2;
            else if(strcmp(atom->name, " C2 ") == 0)
                return C3H1;
            else if(strcmp(atom->name, " N1 ") == 0)
                return N2H0;
            
            // Start common features of pyr ring
            if(strcmp(atom->name, " N9 ") == 0)
                return N3H0;
            else if(strcmp(atom->name, " C8 ") == 0)
                return C3H1;
            else if(strcmp(atom->name, " N7 ") == 0)
                return N2H0;
            else if(strcmp(atom->name, " C6 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " C5 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " C4 ") == 0)
                return C3H0;
            else if(strcmp(atom->name, " N3 ") == 0)
                return N2H0;
            // End common features of pyr ring
            // Start phosphate/sugar backbone
            if(strcmp(atom->name, " P  ") == 0)
                return P4H0;
            else if(strcmp(atom->name, " OP1") == 0)
                return O1H0;
            else if(strcmp(atom->name, " OP2") == 0)
                return O1H0;
            else if(strcmp(atom->name, " O5'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " C5'") == 0)
                return C4H2;
            else if(strcmp(atom->name, " C4'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " C3'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " C2'") == 0)
                if(atom->res_type <= DU)
                    return C4H2;
                else
                    return C4H1;
            else if(strcmp(atom->name, " C1'") == 0)
                return C4H1;
            else if(strcmp(atom->name, " O4'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " O3'") == 0)
                return O2H0;
            else if(strcmp(atom->name, " O2'") == 0)
                return O2H1;
            break;
            // End phosphate/sugar backbone
    }
    
    // TODO: Error?
    return UnknownBonds;
}

void assign_organic_atom_types(Atom* atoms, size_t num_atoms)
{
    for(size_t a = 0; a < num_atoms; a++)
        if(atoms[a].element == Carbon || atoms[a].element == Nitrogen || atoms[a].element == Oxygen ||
                atoms[a].element == Sulfur || atoms[a].element == Phosphorus)
            atoms[a].organic_type = determine_bonding_config(atoms+a);
        else
            atoms[a].organic_type = UnknownBonds;
}

double get_vdw_radius(Atom* atom, bool implicit_hydrogen)
{
    assert(atom->element != UnknownAtom);
    
    if(implicit_hydrogen && atom->organic_type < UnknownBonds)
        return OrganicAtomVdw[atom->organic_type];
    else
        return VdwRadii[atom->element];
}

double get_atom_volume(Atom *atom, bool implicit_hydrogen)
{
    assert(atom->element != UnknownAtom);
    
    if(implicit_hydrogen && atom->organic_type < UnknownBonds)
        return OrganicAtomVolumes[atom->organic_type];
    else
        return VdwSphereVolumes[atom->element];
}

double shrake_rupley_sasa(double* restrict sasa_per_atom, const double probe_radius, const Atom* restrict atoms, const size_t num_atoms,
                          bool implicit_hydrogen, NeighborList* neighborList, coord* restrict sampling_vectors,
                          const size_t num_vectors)
{
    assert(sasa_per_atom != NULL);
    assert(atoms != NULL);
    assert(probe_radius > 0.0);
  
    double total_sasa = 0.0;
    
    // Generate sampling vectors if not provided
    bool local_sampvec = false;
    if(sampling_vectors == NULL)
    {
        // It's either assert this or make num_vectors not const. I'd prefer to ward against bugs for the time being.
        assert(num_vectors % 2 == 1);
        sampling_vectors = malloc(sizeof(coord) * num_vectors);
        golden_spiral_distribution(sampling_vectors, num_vectors);
        local_sampvec = true;
    }
    
    // Generate NeighborList if not provided
    bool local_nl = false;
    if(neighborList == NULL)
    {
        neighborList = malloc(sizeof(NeighborList));
        allocate_make_neighbor_list(neighborList, atoms, num_atoms, probe_radius * 2 + 4.0);
        local_nl = true;
    }
    
    // An array to store the shifted sampling points (gets re-written every loop)
    coord* sampling_points = malloc(sizeof(coord) * num_vectors);
    
    for(size_t i = 0; i < neighborList->cell_lists_len; i++)
    {
        // cells adjacent to `i` are any of the 26 cells with (X, Y, Z) that differ by no more than 1
        // here the first slot is reserved for the self
        size_t adjacent_cells[27] = { i };
        size_t adjacent_count = neighborlist_get_adjacent_cell_indices(adjacent_cells + 1, i, neighborList) + 1;
        
        // for every atom in the current cell...
        for(size_t cell_idx = 0; cell_idx < neighborList->cells_len[i]; cell_idx++)
        {
            size_t atom_idx = neighborList->cell_lists[i][cell_idx];
            Atom atom = atoms[atom_idx];
            // for solvent accessible surface, the test is "What part of a sphere of radius (vdw + probe) is outside of
            // all neighboring atoms sphere's of (radius + probe)?"
            double atom_radius = get_vdw_radius(&atom, implicit_hydrogen);
            double atom_radius_with_probe = atom_radius + probe_radius;
            
            // The count of how many are exposed
            size_t exposed_count = num_vectors;
            
            // Generate sampling points
            for(size_t v = 0; v < num_vectors; v++)
                sampling_points[v] = vec_add(atom.position, vec_multiply(sampling_vectors[v], atom_radius_with_probe));
            
            // for every sampling point on that atom...
            for(size_t v = 0; v < num_vectors; v++)
            {
                coord surface_point = sampling_points[v];
                bool exposed = true;
                // for every cell within range...
                for(size_t j = 0; j < adjacent_count && exposed ; j++)
                {
                    size_t current_cell_idx = adjacent_cells[j];
                    assert(current_cell_idx < neighborList->cell_lists_len);
                    size_t* cell = neighborList->cell_lists[current_cell_idx];
                    // for every neighboring atom in neighboring cells...
                    for(size_t k = 0; k <  neighborList->cells_len[current_cell_idx] && exposed; k++)
                    {
                        size_t neighbor_idx = cell[k];
                        Atom neighbor = atoms[neighbor_idx];
                        double neighbor_radius = get_vdw_radius(&neighbor, implicit_hydrogen) + probe_radius;
                        if(distance_squared(surface_point, neighbor.position) < neighbor_radius * neighbor_radius)
                        {
                            assert(exposed_count > 0);
                            exposed_count--;
                            exposed = false;
                        }
                    }
                }
            }
            
            // SASA is approximately (number of exposed points / total points) * atom surface area
            double sasa_fraction = ((double) exposed_count) / num_vectors;
            // I think this is not the intuitive value to store here, but it is the most useful for SAXS.
            // TODO: Think about how to return the most relevant information here
            sasa_per_atom[atom_idx] = sasa_fraction;
            total_sasa += 4.0 * coord_pi * atom_radius * atom_radius * sasa_fraction;
        }
    }
    
    // Cleanup
    free(sampling_points);
    if(local_sampvec)
        free(sampling_vectors);
    if(local_nl)
    {
        free_neighborlist(neighborList);
        free(neighborList);
    }
    
    return total_sasa;
}
