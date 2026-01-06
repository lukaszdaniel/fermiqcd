/////////////////////////////////////////////////////////////////
/// @file mdp_endianess_converter.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of function swicth_endianess_byte4()
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_ENDIANESS_CONVERTER_
#define MDP_ENDIANESS_CONVERTER_

#include "mdp_macros.h"
#include "mdp_communicator.h"

namespace MDP
{
  /// Converts endianess of object passed by reference
  template <class T>
  void switch_endianess(T &a)
  {
    constexpr size_t size = sizeof(T);

    // Supported size is either 4 or 8 bytes
    if (size == 4 || size == 8)
    {
      char *p = reinterpret_cast<char *>(&a);
      char q[size];

      for (size_t i = 0; i < size; ++i)
      {
        q[i] = p[i];
      }

      for (size_t i = 0; i < size; ++i)
      {
        p[i] = q[size - 1 - i];
      }
    }
    else
    {
      error("switch_endianess: sizeof(T) mismatch");
    }
  }
} // namespace MDP

#endif /* MDP_ENDIANESS_CONVERTER_ */
