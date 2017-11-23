/////////////////////////////////////////////////////////////////
/// @file mdp_global_vars.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// MDP global variables
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef mdp_global_vars_
#define mdp_global_vars_

using namespace std;

typedef unsigned int uint;

//#define CHECK_ALL
//#define MDP_MPI
//#define INCLUDE_DEPRECATED_IO
#define USE_DOUBLE_PRECISION
//#define PARALLEL
//#define NO_POSIX
//#define HAVE_NO_TIMEZONE
//#define SSE2
//#define NO_SSE2_LINALG
//#define DO_NOT_USE_MDP_COMPLEX //define if you want to use standard complex.h header
//#define MDP_NO_LG //define if you want a temporary file to store local-to-global lattice mappings
//#define MATRIXOPTIMIZE
//#define MATRIX_SSE2
//#define AIX
//#define BLOCKSITE 100
//#define TWISTED_BOUNDARY

typedef int mdp_int;

#ifdef USE_DOUBLE_PRECISION
typedef double mdp_real;
#else
typedef float mdp_real;
#endif

const int EVEN = 0;
const int ODD = 1;
const int EVENODD = 2;
const int _NprocMax_ = 256;
double PRECISION = 3.0e-6;
const mdp_int NOWHERE = INT_MAX;

/// Each program should have a name
char *mdp_program_name = (char*) "A generic test program";

/// Filename to store the random seed
char *mdp_random_seed_filename = 0;

/// Used to determine the local endianess of this machine
const unsigned int mdp_local_endianess = 0x87654321;

const double Pi = 3.1415926535897932384626433832795028841971;

/// Set mdp_shutup=true to suppress default output from any part of
/// the program
bool mdp_shutup = true;

/// Default precision used by iterative algorithms such as 
/// mdp_matrix::sin(), mdp_matrix::cos() and mdp_matrix::exp()
double mdp_precision = 1e-5;

void _mpi_error_message(string, string, int);

#endif /* mdp_global_vars_ */
