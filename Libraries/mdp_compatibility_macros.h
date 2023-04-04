/////////////////////////////////////////////////////////////////
/// @file mdp_compatibility_macros.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains macros for backward compatibility now deprecated
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_COMPATIBILITY_MACROS_
#define MDP_COMPATIBILITY_MACROS_

#define myreal mdp_real
#define site mdp_site
#define Complex mdp_complex
#define Matrix mdp_matrix
#define Random mdp_random

// deprecated macros
#define Measure mdp_measure
#define DynamicArray mdp_array
#define JackBoot mdp_jackboot
#define generic_lattice mdp_lattice
#define generic_field mdp_field
#define Matrix_field mdp_matrix_field
#define Vector_field mdp_vector_field
#define NMatrix_field mdp_nmatrix_field
#define mdp_random_generator mdp_prng

#endif /* MDP_COMPATIBILITY_MACROS_ */
