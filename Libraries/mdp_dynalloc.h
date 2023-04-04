/////////////////////////////////////////////////////////////////
/// @file mdp_dynalloc.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Declaration of overloaded new and delete operators
/// to use memalign when compiled with #define SSE2
/// Required for SSE/SSE2 assembly macros
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_DYNALLOC_
#define MDP_DYNALLOC_

#if defined(SSE2) && !defined(OSX)
#include "malloc.h"
void *operator new(size_t size)
{
  void *p = memalign(64, size);
  return p;
}

void operator delete(void *pointer)
{
  free(pointer);
}

void *operator new[](size_t size)
{
  void *p = memalign(64, size);
  return p;
}

void operator delete[](void *pointer)
{
  free(pointer);
}
#endif

#endif /* MDP_DYNALLOC_ */
