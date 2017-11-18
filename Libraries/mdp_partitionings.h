/////////////////////////////////////////////////////////////////
/// @file mdp_partitionings.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Example functions to do parallel partitioning of a lattice
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef mdp_partitionings_
#define mdp_partitionings_

using namespace std;

int default_partitioning0(int *x, int ndim, int *nx) {
	float tpp1 = (float) nx[0] / Nproc;
	int tpp2 = (int) tpp1;
	if (tpp1 - tpp2 > 0)
		tpp2 += 1;
	return x[0] / tpp2;
}

template<int dim>
int default_partitioning(int *x, int ndim, int *nx) {
	float tpp1 = (float) nx[dim] / Nproc;
	int tpp2 = (int) tpp1;
	if (tpp1 - tpp2 > 0)
		tpp2 += 1;
	return x[dim] / tpp2;
}

#endif /* mdp_partitionings_ */
