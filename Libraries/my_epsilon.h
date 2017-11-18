///////////////////////////////////////////////////////////////
/// Full antisymmetric symbol
/// epsilon(i,j,k)
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
