///////////////////////////////////////////////////////////////
/// Full antisymmetric symbol
/// epsilon(i,j,k)

int epsilon(const int& i, const int& j, const int& k) {

	short int sign = -1;

	if (i == k or i == j or j == k)
		sign = 0;
	else if ((i < j and j < k) or (k < i and i < j) or (j < k and k < i))
		sign = 1;

	return sign;
}
////////////////////////////////////////////////////////////////
