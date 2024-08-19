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
  void switch_endianess_byte4(T &a)
  {
    constexpr unsigned int size = 4;
    if (sizeof(T) != size)
      error("switch_endianess_byte: sizeof(T) mismatch");
    char *p = (char *)&a;
    static char q[size];
    for (unsigned int i = 0; i < size; ++i)
    {
      q[i] = p[i];
    }
    for (unsigned int i = 0; i < size; ++i)
    {
      p[i] = q[size - 1 - i];
    }
  }

  template <class T>
  void switch_endianess_byte8(T &a)
  {
    constexpr unsigned int size = 8;
    if (sizeof(T) != size)
      error("switch_endianess_byte: sizeof(T) mismatch");
    char *p = (char *)&a;
    static char q[size];
    for (unsigned int i = 0; i < size; ++i)
    {
      q[i] = p[i];
    }
    for (unsigned int i = 0; i < size; ++i)
    {
      p[i] = q[size - 1 - i];
    }
  }
} // namespace MDP

#endif /* MDP_ENDIANESS_CONVERTER_ */
