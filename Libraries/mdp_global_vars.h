/////////////////////////////////////////////////////////////////
/// @file mdp_global_vars.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// MDP global variables
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_GLOBAL_VARS_
#define MDP_GLOBAL_VARS_

#include <string>
#include <limits>  // for numeric_limits()
#include <cstdint> // for uint8_t

namespace MDP
{
// #define USE_SINGLE_PRECISION // define if you want to use single precision floating point numbers (float) instead of double precision (double)
// #define CHECK_ALL // define if you want to have extra sanity checks (useful for debugging, but slows down the code)
// #define CHECK_BOUNDARY // define if you want to have extra sanity checks for field boundaries (useful for debugging, but slows down the code)
// #define PARALLEL
// #define NO_POSIX
// #define HAVE_NO_TIMEZONE
// #define MDP_USE_SFMT // define if you want to use the SIMD-oriented Fast Mersenne Twister (SFMT) as the random number generator
#define DO_NOT_USE_MDP_COMPLEX // define if you want to use standard complex header
// #define MDP_NO_LG //define if you want a temporary file to store local-to-global lattice mappings
// #define TWISTED_BOUNDARY
// #define LATTICE_DEBUG

#ifdef TWISTED_BOUNDARY
#ifndef BLOCKSITE
#define BLOCKSITE 100
#endif
#endif

    /** @brief Integer value
     *
     * Suitable for large integer values
     *
     * Used for ranges (roughly) +/- 2 bln
     */
    using mdp_int = int32_t;

    /** @brief Short int
     *
     * Suitable for small values
     *
     * Used for ranges -32767, to 32767
     */
    using mdp_sint = int16_t;

    /** @brief Short short int
     *
     * Suitable for small values
     *
     * Used for ranges -128 to 127
     */
    using mdp_ssint = int8_t;

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
     * Used for ranges 0 to 65535
     */
    using mdp_suint = uint16_t;

    /** @brief Short short unsigned int
     *
     * Suitable for small non-negative values
     *
     * Used for ranges 0 to 255
     */
    using mdp_ssuint = uint8_t;

#ifdef USE_SINGLE_PRECISION
    using mdp_real = float;
#else
    using mdp_real = double;
#endif

    enum mdp_parity
    {
        EVEN = 0,
        ODD = 1,
        EVENODD = 2
    };

    constexpr mdp_uint NOWHERE = std::numeric_limits<mdp_uint>::max();

    /// Each program should have a name
    const char *mdp_program_name = "A generic test program";

    /// Filename to store the random seed
    const char *mdp_random_seed_filename = nullptr;

    /// Used to determine the local endianess of this machine
    constexpr mdp_uint mdp_local_endianess = 0x87654321;

    constexpr double Pi = 3.1415926535897932384626433832795028841971;

    /// Default precision used by iterative algorithms such as
    /// mdp_matrix::sin(), mdp_matrix::cos() and mdp_matrix::exp()
    constexpr mdp_real mdp_precision = 1.0 / (1 << 18); // ~3.0e-6 in binary

    void _mpi_error_message(const std::string &message, const std::string &file, int line);
} // namespace MDP

#endif /* MDP_GLOBAL_VARS_ */
