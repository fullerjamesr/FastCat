# FastCat
A tool for high-performance calculation of theoretical small-angle x-ray scattering (SAXS) profiles from input PDB files, and fitting that theoretical profile to experimental data.

## Building and usage
fastcat is written in C11 with no external dependencies outside the C standard library, and so should be straightforward to build. A very basic CMakeLists.txt is provided, but is probably overkill. I develop and test the code against GCC 6.3.0 and ICC 2017 using -std=c11, but I've occasionally used it with GCC 4.9 (which only supports -std=c99) without incident. However, I cannot promise C99 support going forward.

fastcat supports multi-threading via OpenMP, so I recommend enabling that during building if possible.

Usage follows the normal pattern of 
```
$ fastcat [options] structure.pdb [experimental_data.dat]
```
see `fastcat --help` for command line argument information.

## Motivation
This project was concieved as a way to
* Teach myself the C programming language.
* See if I could improve on existing tools in the biomolecular SAXS field.

I was/am comfortable enough with the the technique and underlying math so that I thought I might be able to contribute something useful. In particular, I wanted to see if I could improve the speed of profile generation in order to make structure determination/refinement against SAXS data more feasible, and improve the accuracy of the calcuation in cases where the tradeoff between speed and accuracy was minimal.

## Method citations
* In order to avoid the the poor performance scaling of the Debeye equation as the number of atoms increases, and the accuracy problems of spherical harmonics when applied to anisotropic molecules, I use the algorithm described in Watson and Curtis, "Rapid and accurate calculation of small-angle scattering profiles using the golden ratio" J Appl Cryst (2013). https://doi.org/10.1107/S002188981301666X
* To calculate atomic form factors for each element and maximize their accuracy even at high scattering angles, I use the equations presented in Muhammad and Lee, "New Empirical Equation for the Atomic Form Factor Function in the Momentum Transfer Range, q = 0–50 Å−1 for the Elements in the Range 1≤ Z ≤30" Plos One (2013). https://doi.org/10.1371/journal.pone.0069608
* In order to estimate and correct for excluded volume, particularly in the case of implicit hydrogens, I use the combined atom radii determined by Tsai et al. "The packing density in proteins: standard radii and volumes" J Mol Bio (1999). https://doi.org/10.1006/jmbi.1999.2829
* For modeling the hydration layer, I take the implicit hydration shell approach described by the FoXS program authors: Schneidman-Duhovny et al., "Accurate SAXS Profile Computation and its Assessment by Contrast Variation Experiments" Biophys J (2013). https://dx.doi.org/10.1016/j.bpj.2013.07.020
