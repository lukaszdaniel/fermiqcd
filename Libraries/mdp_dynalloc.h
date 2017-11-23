/////////////////////////////////////////////////////////////////
/// @file mdp_dynalloc.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Declaration of overloaded new and delete operators 
/// to use memalign when compiled with #define SSE2
/// Required for SSE/SSE2 assembly macros
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef mdp_dynalloc_
#define mdp_dynalloc_

using namespace std;

#if defined(SSE2) && !defined(__APPLE__)
#include "malloc.h"
void* operator new(size_t size) {
	void *p=memalign(64,size);
	return p;
}

void operator delete (void* pointer) {
	free(pointer);
}

void* operator new[] (size_t size) {
	void* p=memalign(64,size);
	return p;
}

void operator delete[] (void* pointer) {
	free(pointer);
}
#endif

#endif /* mdp_dynalloc_ */
