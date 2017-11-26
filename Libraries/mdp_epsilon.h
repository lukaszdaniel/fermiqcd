/// @file mdp_epsilon.h
/// @version 2017-11-23
/// @author Lukasz Daniel <lukasz.daniel@gmail.com>
///
/// Contains declaration of the epsilon function
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef my_epsilon_
#define my_epsilon_

using namespace std;

int epsilon(const int& i, const int& j, const int& k) {

	short int sign = -1;

	if (i == k or i == j or j == k)
		sign = 0;
	else if ((i < j and j < k) or (k < i and i < j) or (j < k and k < i))
		sign = 1;

	return sign;
}

#endif /* my_epsilon_ */
