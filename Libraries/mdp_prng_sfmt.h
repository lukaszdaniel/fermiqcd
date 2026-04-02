/////////////////////////////////////////////////////////////////
/// @file mdp_prng_sfmt.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// SIMD-oriented Fast Mersenne Twister (SFMT) random number generator
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PRNG_SFMT_
#define MDP_PRNG_SFMT_

#include <cassert>
#include <cstddef> // for size_t
#include "mdp_global_vars.h"

#define MSK1 0xdfffffefU
#define MSK2 0xddfecb7fU
#define MSK3 0xbffaffffU
#define MSK4 0xbffffff6U
#define PARITY1 0x00000001U
#define PARITY2 0x00000000U
#define PARITY3 0x00000000U
#define PARITY4 0x13c9e684U

namespace MDP
{
  class mdp_prng_sfmt
  {
  private:
    bool initialized;
    static constexpr unsigned int MEXP = 19937;
    static constexpr unsigned int N = 156;
    static constexpr unsigned int N32 = 624;
    static constexpr unsigned int POS1 = 122;
    static constexpr unsigned int SL1 = 18;
    static constexpr unsigned int SL2 = 1;
    static constexpr unsigned int SR1 = 11;
    static constexpr unsigned int SR2 = 1;
    unsigned int idx;

    struct W128_T
    {
      unsigned int u[4];
    };

    using w128_t = struct W128_T;
    w128_t sfmt[N];
    unsigned int *psfmt32;

    /**
     * @brief Ensures that the internal state satisfies the SFMT period property.
     *
     * Performs the period certification step defined by the SFMT algorithm.
     * If the state does not satisfy the required parity condition, the state
     * is modified so that the generator achieves the intended period.
     */
    void period_certification(void)
    {
      static unsigned int parity[4] = {PARITY1, PARITY2, PARITY3, PARITY4};
      unsigned int inner = 0;
      unsigned int work;

      for (unsigned int i = 0; i < 4; i++)
        inner ^= psfmt32[i] & parity[i];
      for (unsigned int i = 16; i > 0; i >>= 1)
        inner ^= inner >> i;
      inner &= 1;
      /* check OK */
      if (inner == 1)
      {
        return;
      }
      /* check NG, and modification */
      for (unsigned int i = 0; i < 4; i++)
      {
        work = 1;
        for (unsigned int j = 0; j < 32; j++)
        {
          if ((work & parity[i]) != 0)
          {
            psfmt32[i] ^= work;
            return;
          }
          work = work << 1;
        }
      }
    }

    /**
     * @brief Performs a right shift of a 128-bit value by a given number of bytes.
     *
     * The input 128-bit value is interpreted as two 64-bit halves and shifted
     * right by (shift * 8) bits. The result is stored in @p out.
     *
     * @param out Pointer to the output 128-bit value.
     * @param in Pointer to the input 128-bit value.
     * @param shift Number of bytes to shift.
     */
    void rshift128(w128_t *out, w128_t const *in, int shift)
    {
      size_t th = ((size_t)in->u[3] << 32) | ((unsigned int)in->u[2]);
      size_t tl = ((size_t)in->u[1] << 32) | ((unsigned int)in->u[0]);

      size_t oh = th >> (shift * 8);
      size_t ol = tl >> (shift * 8);
      ol |= th << (64 - shift * 8);

      out->u[1] = (unsigned int)(ol >> 32);
      out->u[0] = (unsigned int)ol;
      out->u[3] = (unsigned int)(oh >> 32);
      out->u[2] = (unsigned int)oh;
    }

    /**
     * @brief Performs a left shift of a 128-bit value by a given number of bytes.
     *
     * The input 128-bit value is interpreted as two 64-bit halves and shifted
     * left by (shift * 8) bits. The result is stored in @p out.
     *
     * @param out Pointer to the output 128-bit value.
     * @param in Pointer to the input 128-bit value.
     * @param shift Number of bytes to shift.
     */
    void lshift128(w128_t *out, w128_t const *in, int shift)
    {
      size_t th = ((size_t)in->u[3] << 32) | ((unsigned int)in->u[2]);
      size_t tl = ((size_t)in->u[1] << 32) | ((unsigned int)in->u[0]);

      size_t oh = th << (shift * 8);
      size_t ol = tl << (shift * 8);
      oh |= tl >> (64 - shift * 8);

      out->u[1] = (unsigned int)(ol >> 32);
      out->u[0] = (unsigned int)ol;
      out->u[3] = (unsigned int)(oh >> 32);
      out->u[2] = (unsigned int)oh;
    }

    /**
     * @brief Generates a new block of pseudorandom numbers for the entire state.
     *
     * Updates the internal SFMT state array using the recursion function.
     * This function is called when the current state buffer has been exhausted.
     */
    void gen_rand_all(void)
    {
      mdp_uint i;
      w128_t *r1 = &sfmt[N - 2];
      w128_t *r2 = &sfmt[N - 1];

      for (i = 0; i < N - POS1; i++)
      {
        do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1], r1, r2);
        r1 = r2;
        r2 = &sfmt[i];
      }
      for (; i < N; i++)
      {
        do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1 - N], r1, r2);
        r1 = r2;
        r2 = &sfmt[i];
      }
    }

    /**
     * @brief Core SFMT recursion step.
     *
     * Computes a new 128-bit state value using the SFMT recurrence relation
     * based on four previous state values.
     *
     * @param r Output state value.
     * @param a Input state value A.
     * @param b Input state value B.
     * @param c Input state value C.
     * @param d Input state value D.
     */
    void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *c,
                      w128_t *d)
    {
      w128_t x;
      w128_t y;

      lshift128(&x, a, SL2);
      rshift128(&y, c, SR2);

      r->u[0] = a->u[0] ^ x.u[0] ^ ((b->u[0] >> SR1) & MSK1) ^ y.u[0] ^ (d->u[0] << SL1);
      r->u[1] = a->u[1] ^ x.u[1] ^ ((b->u[1] >> SR1) & MSK2) ^ y.u[1] ^ (d->u[1] << SL1);
      r->u[2] = a->u[2] ^ x.u[2] ^ ((b->u[2] >> SR1) & MSK3) ^ y.u[2] ^ (d->u[2] << SL1);
      r->u[3] = a->u[3] ^ x.u[3] ^ ((b->u[3] >> SR1) & MSK4) ^ y.u[3] ^ (d->u[3] << SL1);
    }

    /**
     * @brief Returns the next 32-bit pseudorandom integer.
     *
     * If the internal buffer is exhausted, a new batch of random values
     * is generated before returning the next value.
     *
     * @return 32-bit pseudorandom unsigned integer.
     */
    unsigned int gen_rand32()
    {
      unsigned int r;
      assert(initialized);

      if (idx >= N32)
      {
        gen_rand_all();
        idx = 0;
      }

      r = psfmt32[idx++];
      return r;
    }

  public:
    mdp_prng_sfmt(mdp_int k = 0)
    {
      if (k == 0)
        initialize(ME);
    }

    /**
     * @brief Initializes the generator with a given seed.
     *
     * Fills the internal state using the standard SFMT initialization
     * procedure and performs period certification.
     *
     * @param seed Initial seed value.
     */
    void initialize(unsigned int seed)
    {
      psfmt32 = (unsigned int *)&(sfmt[0].u[0]);

      psfmt32[0] = seed;
      for (unsigned int i = 1; i < N32; i++)
      {
        psfmt32[i] = 1812433253UL * (psfmt32[i - 1] ^ (psfmt32[i - 1] >> 30)) + i;
      }

      idx = N32;
      period_certification();
      initialized = true;
    }

    /**
     * @brief Returns a pseudorandom floating-point value in the range [0, 1].
     *
     * The value is generated from a 32-bit random integer and normalized
     * to the interval [0, 1].
     *
     * @return Random float in the range [0, 1].
     */
    inline float plain() noexcept
    {
      unsigned int v = gen_rand32();
      return v * (1.0 / 4294967295.0);
    }
  };
} // namespace MDP

#endif /* MDP_PRNG_SFMT_ */
