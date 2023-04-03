/////////////////////////////////////////////////////////////////
/// @file mdp.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Includes all mdp_*.h header files
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef MDP_
#define MDP_

// ///////////////////////////////////////////////////////////////////////////
// include the usual libraries (works on gcc and VC)
// ///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cassert>
#include <typeinfo>
#include <cstdint>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <climits>
#ifndef _WIN64
#include "glob.h"
#endif
#ifndef NO_POSIX
#include <unistd.h>
#include <sys/time.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN64
#include "sys/socket.h"
#endif
#include <fcntl.h>
#endif

// ///////////////////////////////////////////////////////////////////////////
// this file includes the version number
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_version.h"

// ///////////////////////////////////////////////////////////////////////////
// all global macros used by MDP
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_macros.h"

// ///////////////////////////////////////////////////////////////////////////
// all global varibales except mdp,mpi and mdp_random
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_global_vars.h"

// ///////////////////////////////////////////////////////////////////////////
// faster dynamic allocation (no exceptions)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_dynalloc.h"

// ///////////////////////////////////////////////////////////////////////////
// function to convert endianess
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_endianess_converter.h"

// ///////////////////////////////////////////////////////////////////////////
// this is the official mdp_timer (replaces JIM_timer since not portable
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_timer.h"

// ///////////////////////////////////////////////////////////////////////////
// mdp implementation of complex numbers (portable ansi)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_complex.h"

// ///////////////////////////////////////////////////////////////////////////
// delta function
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_delta.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of multidimentional array, better than STL
// (class mdp_array used to be class DynamicArray)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_array.h"

// ///////////////////////////////////////////////////////////////////////////
// this file includes stuff for sse2 optimization
// ///////////////////////////////////////////////////////////////////////////
#ifdef SSE2
#include "fermiqcd_sse.h"
#endif

// ///////////////////////////////////////////////////////////////////////////
// implementation of the mdp_matrix object
// (class mdp_matrix is a more general implementation of class mdp_matrix)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_matrix.h"

// ///////////////////////////////////////////////////////////////////////////
// class for logging (not very developed!)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_log.h"

// ///////////////////////////////////////////////////////////////////////////
// parallel simulator
// handy for debugging, multithreading and mosix systems
// ///////////////////////////////////////////////////////////////////////////
#ifndef NO_POSIX
#include "mdp_psim.h"
#endif

// ///////////////////////////////////////////////////////////////////////////
// this is a wrapper to Message Passing Interface (is one uses it)
// replace this functions to change communication protocol
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_communicator.h"

// ///////////////////////////////////////////////////////////////////////////
// this defined the class mdp_prng and the obj mdp_random
// (attention that ::SU<T>(int n) only works with gcc,
//  VC does not support templates, therefore ::SU(int n) only for float)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_prng.h"

// ///////////////////////////////////////////////////////////////////////////
// mdp_jackboot is a class for statistical analysis
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_jackboot.h"

// ///////////////////////////////////////////////////////////////////////////
// this is a collection of possible lattice topologies
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_topologies.h"

// ///////////////////////////////////////////////////////////////////////////
// a collection of possible lattice partitioning
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_partitionings.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the class mdp_lattice
// (used to be generic_lattice)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_lattice.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of a vector on a lattice (used internally for conversions)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_vector.h"

// ///////////////////////////////////////////////////////////////////////////
// class mdp_site (used to be site) a point on a lattice
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_site.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the class mdp_field<>
// (used to be generic_field<>)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_field.h"

// ///////////////////////////////////////////////////////////////////////////
// various other utilities
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_utils.h"
#include "mdp_postscript.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the communicaton function mdp_field::update()
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_field_update.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the I/O function mdp_field::load()
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_field_load.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the I/O function mdp_field::save()
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_field_save.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the I/O function mdp_field::save_vtk()
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_field_save_vtk.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the I/O function mdp_field::save_vtk()
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_save_partitioning_vtk.h"

#ifdef INCLUDE_DEPRECATED_IO
#include "mdp_deprecatedIO.h"
#endif

// ///////////////////////////////////////////////////////////////////////////
// an auxiliary function that returns (-1)^n
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_mod2sign.h"

// ///////////////////////////////////////////////////////////////////////////
// very clever function to compute permutations of lists
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_permutations.h"

// ///////////////////////////////////////////////////////////////////////////
// implementation of the class mdp_complex_field
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_complex_field.h"

// ///////////////////////////////////////////////////////////////////////////
// an mdp_field of matrices
// (more general than mdp_matrix_field)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_matrix_field.h"

// ///////////////////////////////////////////////////////////////////////////
// an mdp_field of a vector matrices
// (more general than Nmdp_matrix_field)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_vector_field.h"

// ///////////////////////////////////////////////////////////////////////////
// an mdp_field of vector (1xN matrix)
// (more general than Vector_field)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_nmatrix_field.h"

// ///////////////////////////////////////////////////////////////////////////
// compatibility functions map MDP 1.3 into MDP 2.0 or higher
// (only the syntax of mdp_random::SU(int n) is not portable
//  everything else is portable if according to specs)
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_compatibility_macros.h"

// ///////////////////////////////////////////////////////////////////////////
// functions to prompt the user for input values
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_prompt.h"

// ///////////////////////////////////////////////////////////////////////////
// a container class for experimental results
// ///////////////////////////////////////////////////////////////////////////
#include "mdp_measure.h"

#include "mdp_fitting_functions.h"
#include "mdp_header.h"
#include "mdp_nvector_field.h"
#include "mdp_prng_sfmt.h"
#include "mdp_swap.h"

#endif /* MDP_ */