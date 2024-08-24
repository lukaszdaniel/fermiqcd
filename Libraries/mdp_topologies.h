/////////////////////////////////////////////////////////////////
/// @file mdp_topologies.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Examples of lattice topologies
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_TOPOLOGIES_
#define MDP_TOPOLOGIES_

namespace MDP
{
	// //////////////////////////////////////////////////////
	// Basic topologies:
	// //////////////////////////////////////////////////////

	void torus_topology(const int mu,
						int x_dw[],
						const int x[],
						int x_up[],
						const int ndim,
						const int nx[])
	{
		for (int nu = 0; nu < ndim; nu++)
			if (nu == mu)
			{
				x_dw[mu] = (x[mu] - 1 + nx[mu]) % nx[mu];
				x_up[mu] = (x[mu] + 1) % nx[mu];
			}
			else
				x_up[nu] = x_dw[nu] = x[nu];
	}

	void box_topology(const int mu,
					  int x_dw[],
					  const int x[],
					  int x_up[],
					  const int ndim,
					  const int nx[])
	{
		torus_topology(mu, x_dw, x, x_up, ndim, nx);
		if (x[mu] == 0)
			x_dw[mu] = x[mu];
		if (x[mu] == nx[mu] - 1)
			x_up[mu] = x[mu];
	}

	void moebious_topology(const int mu,
						   int x_dw[],
						   const int x[],
						   int x_up[],
						   const int ndim,
						   const int nx[])
	{
		torus_topology(mu, x_dw, x, x_up, ndim, nx);
		if (mu == 0)
		{
			if (x[0] == 0)
				x_dw[1] = nx[1] - x[1] - 1;
			if (x[0] == nx[0] - 1)
				x_up[1] = nx[1] - x[1] - 1;
		}
	}

	void open_cylinder(const int mu,
					   int x_dw[],
					   const int x[],
					   int x_up[],
					   const int ndim,
					   const int nx[])
	{
		torus_topology(mu, x_dw, x, x_up, ndim, nx);
		if ((mu == 0) && (x[0] == 0))
			x_dw[0] = x[0];
		if ((mu == 0) && (x[0] == nx[0] - 1))
			x_up[0] = x[0];
	}
} // namespace MDP

#endif /* MDP_TOPOLOGIES_ */
