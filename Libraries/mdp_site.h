/////////////////////////////////////////////////////////////////
/// @file mdp_site.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains delcaration of class mdp_site for complex numbers
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_SITE_
#define MDP_SITE_

namespace MDP
{
  /// checks which process of the lattice a stores locally the site of
  /// coordinates x0,x1,x2,...,x9
  /// to be used before calling mdp_site::set()
  /// (note: prototyping of friend functions is required by some compilers)
  int on_which_process(mdp_lattice &a,
                       int x0 = 0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
                       int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0);

  /// @brief site object to loop on a lattice
  ///
  /// Example:
  /// @verbatim
  ///   int box[]={10,10,10};
  ///   mdp_lattice lattice(3,box);
  ///   mdp_site x(lattice);
  ///   forallsites(x) cout << x << endl;
  ///   if(on_which_process(lattice,1,1,1)==ME) {
  ///      x.set(1,1,1);
  ///      cout << lattice.random(x).plain() << endl;
  ///   }
  /// @endverbatim
  class mdp_site
  {
  private:
    mdp_lattice *ptr; /// this points to the lattice for this field

    mdp_site()
    {
      // one whould not use this!
      ptr = nullptr;
      idx = 0;
    }

    void on(const mdp_lattice &a)
    {
      ptr = (mdp_lattice *)&a;
      idx = (*ptr).start[ME][0];
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = 0;
#endif
    }

  public:
    mdp_int idx; /// value of the mdp_site in local coordinate
#ifdef BLOCKSITE
    int block[BLOCKSITE];
#endif

    /// declares object of class mdp_site living on the
    /// lattice passed by reference
    mdp_site(const mdp_lattice &a)
    {
      on(a);
    }

    /// returns by reference the lattice the site lives on
    mdp_lattice &lattice()
    {
      return *ptr;
    }

    mdp_site(mdp_int i, mdp_lattice *ptr2)
    {
      idx = i;
      ptr = ptr2;
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = 0;
#endif
    }

#ifdef BLOCKSITE
    mdp_site(mdp_int i, mdp_lattice *ptr2, int b[], int sign = 0, int mu = 0)
    {
      idx = i;
      ptr = ptr2;
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = b[k];
      block[mu] += sign;
    }
#endif
    mdp_site(const mdp_site &x)
    {
      idx = x.idx;
      ptr = x.ptr;
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = x.block[k];
#endif
    }

    mdp_site operator=(mdp_int i)
    {
      idx = lattice().start[ME][0] + i;
      return mdp_site(idx, ptr);
    }

    mdp_site operator=(mdp_site x)
    {
      if (ptr == x.ptr)
        idx = x.idx;
      else
        set_global(x.global_index());
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = x.block[k];
      return mdp_site(idx, ptr, block);
#endif
      return mdp_site(idx, ptr);
    }

    bool operator==(mdp_site x)
    {
      if ((idx == NOWHERE) || (x.idx == NOWHERE))
        return false;
      return (global_index() == x.global_index());
    }

    bool operator!=(mdp_site x)
    {
      return !(*this == x);
    }

    void start(int np = 0)
    {
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        block[k] = 0;
#endif
      idx = lattice().start[ME][np];
    }

    void next()
    {
      ++idx;
    }

    /** @brief checks if the site is inside the portion of the lattice stored by
     * the current process
     */
    bool is_in()
    {
      return ((idx >= lattice().start[ME][0]) && (idx < lattice().stop[ME][1]));
    }

    /** @brief checks if the site is inside the portion of the lattice stored by
     * the current process or if the site is in a local copy of a remote
     * site
     */
    int is_here()
    {
      return ((idx >= 0) && (idx < lattice().enclosing_volume()));
    }

    /** @brief returns the parity EVEN or ODD of the site
     */
    int parity()
    {
      return lattice().site_parity(idx);
    }

    /** @brief true if the site is stored locally as a copy of a site local
     * in another process
     */
    int is_in_boundary()
    {
      return (lattice().wh[idx] != ME);
    }

    /// returns the local index of the site
    /// local index is assigned by the process to the local sites and copies
    /// of remote sites. local index is not unique thoughout the lattice.
    mdp_int local_index()
    {
      return idx;
    }

    /** @brief returns the global (unique) index of the site
     */
    mdp_int global_index()
    {
      return lattice().global(idx);
    }

    /** @brief sets the site by its local index (dangerous)
     */
    void set_local(mdp_int idx2)
    {
      idx = idx2;
    }

    /** @brief sets the site by its global index
     */
    void set_global(mdp_int idx_gl)
    {
      mdp_int idx2 = lattice().local(idx_gl);
      if (idx2 == NOWHERE)
        error("set_global() trying to access a site that is not here");
      idx = idx2;
    }

    /** @brief returns the site shifted forward in direction mu=(0...ndim-1)
     */
    mdp_site operator+(int mu)
    {
      mdp_int idx2 = lattice().up[idx][mu];
      if (idx2 == NOWHERE)
      {
        std::cout << ME << " " << (*this)(0) << (*this)(1)
                  << (*this)(2) << (*this)(3) << " " << mu << std::endl;
        error("You cannot exit from your portion of the lattice");
      }
#ifdef BLOCKSITE
      if ((mu < BLOCKSITE) && (lattice().co[idx][mu] == lattice().nx[mu] - 1))
        return mdp_site(idx2, ptr, block, 1, mu);
      return mdp_site(idx2, ptr, block);
#endif
      return mdp_site(idx2, ptr);
    }

    /** @brief returns the site shifted backwards in direction mu=(0...ndim-1)
     */
    mdp_site operator-(int mu)
    {
      mdp_int idx2 = lattice().dw[idx][mu];
      if (idx2 == NOWHERE)
        error("You cannot exit from your portion of the lattice");
#ifdef BLOCKSITE
      if ((mu < BLOCKSITE) && (lattice().co[idx][mu] == 0))
        return mdp_site(idx2, ptr, block, -1, mu);
      return mdp_site(idx2, ptr, block);
#endif
      return mdp_site(idx2, ptr);
    }

    /** @brief returns a site shifted i position (backwards if i<0 or forward if i>0)
     * in direction mu=(0...mdim-1)
     */
    mdp_site hop(int i, int mu)
    {
      mdp_site y(lattice());
      y = (*this);
      while (i != 0)
      {
        if (i < 0)
        {
          y = y - mu;
          i++;
        }
        if (i > 0)
        {
          y = y + mu;
          i--;
        }
      }
      return y;
    }

    /** @brief sets the site to the coordinates stored in vector v
     */
    mdp_site operator=(mdp_vector v)
    {
      set(v[0], v[1], v[2], v[3], v[4],
          v[5], v[6], v[7], v[8], v[9]);
      return *this;
    }

    /// returns a site similar to the present but
    /// each coordinates mu of the site shifted according to v[mu]
    mdp_site operator+(mdp_vector v)
    {
      int mu, step;
      mdp_site y = *this;
      for (mu = 0; mu < lattice().n_dimensions(); mu++)
      {
        if (v[mu] > 0)
          for (step = 0; step < v[mu]; step++)
            y = y + mu;
        else
          for (step = 0; step < -v[mu]; step++)
            y = y - mu;
      }
      return y;
    }

    /// returns a site similar to the present but
    /// each coordinates mu of the site shifted according to -v[mu]
    mdp_site operator-(mdp_vector v)
    {
      int mu, step;
      mdp_site y = *this;
      for (mu = 0; mu < lattice().n_dimensions(); mu++)
      {
        if (v[mu] > 0)
          for (step = 0; step < v[mu]; step++)
            y = y - mu;
        else
          for (step = 0; step < -v[mu]; step++)
            y = y + mu;
      }
      return y;
    }

    /// returns mu coordinate of the site
    int operator()(int mu)
    {
      return lattice().co[idx][mu];
    }

    void operator=(int *x)
    {
      int ndim = lattice().n_dimensions();
      idx = x[0];
      for (int mu = 1; mu < ndim; mu++)
        idx = idx * lattice().size(mu) + x[mu];
      idx = lattice().local(idx);
      if (idx == NOWHERE)
      {
        bool print = mdp.print;
        mdp.print = true;
        mdp << "Warning message from ME=" << ME << ":\n";
        mdp << "You assigned a site that is not here!\n";
        mdp.print = print;
      }
    }

    /// sets the site to a the location specified by the coordinates
    /// and assumes the site is local (or at least a copy)
    /// @see on_which_process()
    void set(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
             int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
    {
      int ndim = lattice().n_dimensions();
      idx = x0;
      if (ndim > 1)
        idx = idx * lattice().size(1) + x1;
      if (ndim > 2)
        idx = idx * lattice().size(2) + x2;
      if (ndim > 3)
        idx = idx * lattice().size(3) + x3;
      if (ndim > 4)
        idx = idx * lattice().size(4) + x4;
      if (ndim > 5)
        idx = idx * lattice().size(5) + x5;
      if (ndim > 6)
        idx = idx * lattice().size(6) + x6;
      if (ndim > 7)
        idx = idx * lattice().size(7) + x7;
      if (ndim > 8)
        idx = idx * lattice().size(8) + x8;
      if (ndim > 9)
        idx = idx * lattice().size(9) + x9;

      idx = lattice().local(idx);

      if (idx == NOWHERE)
      {
        bool print = mdp.print;
        mdp.print = true;
        mdp << "Warning message from ME=" << ME << ":\n";
        mdp << "You assigned a site that is not here!\n";
        mdp.print = print;
      }
    }

    bool operator==(const int *x)
    {
      int ndim = lattice().n_dimensions();
      bool is_it = true;
      for (int mu = 0; mu < ndim; mu++)
        if (x[mu] != lattice().co[idx][mu])
          is_it = false;
      return is_it;
    }

    bool operator!=(const int *x)
    {
      return !(*this == x);
    }

    /// checks the site coordinates vs the coordinates passed as args
    int is_equal(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
                 int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
    {
      int ndim = lattice().n_dimensions();

      if ((ndim > 0) && (x0 != lattice().co[idx][0]))
        return false;
      if ((ndim > 1) && (x1 != lattice().co[idx][1]))
        return false;
      if ((ndim > 2) && (x2 != lattice().co[idx][2]))
        return false;
      if ((ndim > 3) && (x3 != lattice().co[idx][3]))
        return false;
      if ((ndim > 4) && (x4 != lattice().co[idx][4]))
        return false;
      if ((ndim > 5) && (x5 != lattice().co[idx][5]))
        return false;
      if ((ndim > 6) && (x6 != lattice().co[idx][6]))
        return false;
      if ((ndim > 7) && (x7 != lattice().co[idx][7]))
        return false;
      if ((ndim > 8) && (x8 != lattice().co[idx][8]))
        return false;
      if ((ndim > 9) && (x9 != lattice().co[idx][9]))
        return false;
      return true;
    }

    /// converts a site into a binary number
    /// to be used only if the site is a vertex of an hypercube centered
    /// at the origin. this is used to make staggered mesons
    friend mdp_int site2binary(mdp_site x)
    {
      int a = 0;
      for (int mu = 0; mu < x.lattice().n_dimensions(); mu++)
      {
#ifdef CHECK_ALL
        if (fabs(0.5 - x(mu)) > 1)
          error("site2binary");
#endif
        a += (0x1 << mu) * x(mu);
      }
      return a;
    }

    /// checks which process of the lattice a stores locally the site of
    /// coordinates x0,x1,x2,...,x9
    /// to be used before calling mdp_site::set()
    friend int on_which_process(mdp_lattice &a, int x0, int x1, int x2, int x3,
                                int x4, int x5, int x6, int x7, int x8, int x9)
    {
      int x[10];
      x[0] = x0;
      if (a.n_dimensions() > 1)
        x[1] = x1;
      if (a.n_dimensions() > 2)
        x[2] = x2;
      if (a.n_dimensions() > 3)
        x[3] = x3;
      if (a.n_dimensions() > 4)
        x[4] = x4;
      if (a.n_dimensions() > 5)
        x[5] = x5;
      if (a.n_dimensions() > 6)
        x[6] = x6;
      if (a.n_dimensions() > 7)
        x[7] = x7;
      if (a.n_dimensions() > 8)
        x[8] = x8;
      if (a.n_dimensions() > 9)
        x[9] = x9;
      return (*(a.where))(x, a.n_dimensions(), a.nx);
    }
  };

#ifdef MDP_LATTICE

  /// Returns the local object mdp_prng at site x of the lattice
  mdp_prng &mdp_lattice::random(mdp_site x)
  {
    if (local_random_generator)
    {
      if (!x.is_in())
        error("request the random generator of a non local site");
      return random_obj[x.idx - start[ME][0]];
    }
    return mdp_random;
  }

#endif

  /// When compiled with TWISTED_BOUNDARY the mdp_site class keeps track of
  /// sites that moved around the boundary of the torus topology. this function
  /// Returns false if this is one such site, true otherwise.
  int in_block(mdp_site x)
  {
#ifdef TWISTED_BOUNDARY
    for (int mu = 0; mu < x.lattice().n_dimensions(); mu++)
      if (x.block[mu] != 0)
        return false;
#endif
    return true;
  }

  std::ostream &operator<<(std::ostream &os, mdp_site &x)
  {
    for (int i = 0; i < x.lattice().n_dimensions(); i++)
    {
      if (i == 0)
        os << "(" << x(i);
      else
        os << "," << x(i);
    }
    os << ")";
    return os;
  }
} // namespace MDP

#endif /* MDP_SITE_ */
