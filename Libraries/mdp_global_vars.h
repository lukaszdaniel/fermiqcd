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

typedef unsigned int uint;

const int EVEN = 0;
const int ODD = 1;
const int EVENODD = 2;
const int _NprocMax_ = 256;
double PRECISION = 3.0e-6;

/// Each program should have a name
char *mdp_program_name = (char *)"A generic test program";

/// Filename to store the random seed
char *mdp_random_seed_filename = 0;

/// Used to determine the local endianess of this machine
const unsigned int mdp_local_endianess = 0x87654321;

const double Pi = 3.1415926535897932384626433832795028841971;

/// Set mdp_shutup=true to suppress default output from any part of
/// The program
bool mdp_shutup = false;

/// Default precision used by iterative algorithms such as
/// mdp_matrix::sin(), mdp_matrix::cos() and mdp_matrix::exp()
double mdp_precision = 1e-5;

#ifdef USE_DOUBLE_PRECISION
typedef double mdp_real;
#else
typedef float mdp_real;
#endif

/** @brief Integer value
 *
 * Suitable for large integer values
 *
 * Used for ranges (roughly) +/- 2 bln
 */
using mdp_int = int32_t;

/** @brief unsigned int
 *
 * Suitable for large non-negative values
 *
 * Used for ranges from 0 to (roughly) 4 bln
 */
using mdp_uint = uint32_t;

/** @brief Short unsigned int
 *
 * Suitable for small non-negative values
 *
 * Used for ranges 0 to 255
 */
using mdp_suint = uint8_t;

void _mpi_error_message(std::string, std::string, int);
