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
    mdp_lattice *m_lattice; /// this points to the lattice for this field
    mdp_int m_idx;          /// value of the mdp_site in local coordinate

#ifdef BLOCKSITE
    int m_block[BLOCKSITE];
#endif

  public:
    /** @brief declares object of class mdp_site living on the
     * lattice passed by reference
     */
    mdp_site(const mdp_lattice &a)
    {
      m_lattice = (mdp_lattice *)&a;
      m_idx = (*m_lattice).start0(ME, 0);
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = 0;
#endif
    }

    mdp_site(mdp_int i, mdp_lattice *ptr2)
    {
      m_idx = i;
      m_lattice = ptr2;
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = 0;
#endif
    }

#ifdef BLOCKSITE
    mdp_site(mdp_int i, mdp_lattice *ptr2, int b[], int sign = 0, int mu = 0)
    {
      m_idx = i;
      m_lattice = ptr2;
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = b[k];

      if (mu < BLOCKSITE)
      {
        m_block[mu] += sign;
      }
      else
      {
        error("BLOCKSITE < lattice dimension");
      }
    }
#endif

    /** @brief Copy constructor
     */
    mdp_site(const mdp_site &x)
    {
      m_idx = x.m_idx;
      m_lattice = x.m_lattice;
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = x.m_block[k];
#endif
    }

    /** @brief returns by reference the lattice the site lives on
     */
    mdp_lattice &lattice()
    {
      return *m_lattice;
    }

#ifdef BLOCKSITE
    int block(mdp_int mu) const
    {
      if (mu < BLOCKSITE)
      {
        return m_block[mu];
      }
      else
      {
        error("BLOCKSITE < lattice dimension");
      }
    }
#endif

    mdp_site operator=(mdp_int i)
    {
      m_idx = lattice().start0(ME, 0) + i;
      return mdp_site(m_idx, m_lattice);
    }

    mdp_site operator=(mdp_site x)
    {
      if (m_lattice == x.m_lattice)
      {
        m_idx = x.m_idx;
      }
      else
      {
        set_global(x.global_index());
      }

#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = x.m_block[k];
      return mdp_site(m_idx, m_lattice, m_block);
#endif

      return mdp_site(m_idx, m_lattice);
    }

    bool operator==(mdp_site x)
    {
      if ((m_idx == NOWHERE) || (x.m_idx == NOWHERE))
        return false;
      return (global_index() == x.global_index());
    }

    bool operator!=(mdp_site x)
    {
      return !(*this == x);
    }

    void start(mdp_int np = 0)
    {
#ifdef BLOCKSITE
      for (mdp_int k = 0; k < BLOCKSITE; k++)
        m_block[k] = 0;
#endif
      m_idx = lattice().start0(ME, np);
    }

    void next()
    {
      ++m_idx;
    }

    /** @brief checks if the site is inside the portion of the lattice stored by
     * the current process
     */
    bool is_in()
    {
      return ((m_idx >= lattice().start0(ME, 0)) && (m_idx < lattice().stop0(ME, 1)));
    }

    /** @brief checks if the site is inside the portion of the lattice stored by
     * the current process or if the site is in a local copy of a remote
     * site
     */
    int is_here()
    {
      return ((m_idx >= 0) && (m_idx < lattice().enclosing_volume()));
    }

    /** @brief returns the parity EVEN or ODD of the site
     */
    mdp_int parity()
    {
      return lattice().site_parity(m_idx);
    }

    /** @brief true if the site is stored locally as a copy of a site local
     * in another process
     */
    bool is_in_boundary()
    {
      return (lattice().process(m_idx) != ME);
    }

    /** @brief returns the local index of the site
     *
     * local index is assigned by the process to the local sites and copies
     * of remote sites. local index is not unique thoughout the lattice.
     */
    mdp_int local_index() const
    {
      return m_idx;
    }

    /** @brief returns the global (unique) index of the site
     */
    mdp_int global_index()
    {
      return lattice().global(m_idx);
    }

    /** @brief sets the site by its local index (dangerous)
     */
    void set_local(mdp_int idx2)
    {
      m_idx = idx2;
    }

    /** @brief sets the site by its global index
     */
    void set_global(mdp_int idx_gl)
    {
      mdp_int idx2 = lattice().local(idx_gl);
      if (idx2 == NOWHERE)
        error("set_global() trying to access a site that is not here");
      m_idx = idx2;
    }

    /** @brief returns the site shifted forward in direction mu=(0...ndim-1)
     */
    mdp_site operator+(mdp_int mu)
    {
      mdp_int idx2 = lattice().move_up(m_idx, mu);
      if (idx2 == NOWHERE)
      {
        std::cout << ME << " " << (*this)(0) << (*this)(1)
                  << (*this)(2) << (*this)(3) << " " << mu << std::endl;
        error("You cannot exit from your portion of the lattice");
      }
#ifdef BLOCKSITE
      if ((mu < BLOCKSITE) && (lattice().coordinate(m_idx, mu) == lattice().size(mu) - 1))
      {
        return mdp_site(idx2, m_lattice, m_block, 1, mu);
      }
      return mdp_site(idx2, m_lattice, m_block);
#endif

      return mdp_site(idx2, m_lattice);
    }

    /** @brief returns the site shifted backwards in direction mu=(0...ndim-1)
     */
    mdp_site operator-(int mu)
    {
      mdp_int idx2 = lattice().move_down(m_idx, mu);
      if (idx2 == NOWHERE)
        error("You cannot exit from your portion of the lattice");
#ifdef BLOCKSITE
      if ((mu < BLOCKSITE) && (lattice().coordinate(m_idx, mu) == 0))
      {
        return mdp_site(idx2, m_lattice, m_block, -1, mu);
      }
      return mdp_site(idx2, m_lattice, m_block);
#endif

      return mdp_site(idx2, m_lattice);
    }

    /** @brief returns a site shifted i position (backwards if i<0 or forward if i>0)
     * in direction mu=(0...mdim-1)
     */
    mdp_site hop(int i, mdp_int mu)
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
      set(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]);
      return *this;
    }

    /** @brief returns a site similar to the present but
     * each coordinates mu of the site shifted according to v[mu]
     */
    mdp_site operator+(mdp_vector v)
    {
      mdp_site y = *this;
      for (mdp_int mu = 0; mu < lattice().n_dimensions(); mu++)
      {
        if (v[mu] > 0)
          for (mdp_int step = 0; step < v[mu]; step++)
            y = y + mu;
        else
          for (mdp_int step = 0; step < -v[mu]; step++)
            y = y - mu;
      }
      return y;
    }

    /** @brief returns a site similar to the present but
     * each coordinates mu of the site shifted according to -v[mu]
     */
    mdp_site operator-(mdp_vector v)
    {
      mdp_site y = *this;
      for (mdp_int mu = 0; mu < lattice().n_dimensions(); mu++)
      {
        if (v[mu] > 0)
          for (mdp_int step = 0; step < v[mu]; step++)
            y = y - mu;
        else
          for (mdp_int step = 0; step < -v[mu]; step++)
            y = y + mu;
      }
      return y;
    }

    /// returns mu coordinate of the site
    mdp_int operator()(mdp_int mu)
    {
      return lattice().coordinate(m_idx, mu);
    }

    void operator=(int *x)
    {
      mdp_int ndim = lattice().n_dimensions();
      m_idx = x[0];
      for (mdp_int mu = 1; mu < ndim; mu++)
        m_idx = m_idx * lattice().size(mu) + x[mu];
      m_idx = lattice().local(m_idx);
      if (m_idx == NOWHERE)
      {
        bool print = mdp.printing();
        mdp.enablePrinting();
        mdp << "Warning message from ME=" << ME << ":\n";
        mdp << "You assigned a site that is not here!\n";
        mdp.restorePrinting(print);
      }
    }

    /// sets the site to a the location specified by the coordinates
    /// and assumes the site is local (or at least a copy)
    /// @see on_which_process()
    void set(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
             int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
    {
      mdp_int ndim = lattice().n_dimensions();
      m_idx = x0;
      if (ndim > 1)
        m_idx = m_idx * lattice().size(1) + x1;
      if (ndim > 2)
        m_idx = m_idx * lattice().size(2) + x2;
      if (ndim > 3)
        m_idx = m_idx * lattice().size(3) + x3;
      if (ndim > 4)
        m_idx = m_idx * lattice().size(4) + x4;
      if (ndim > 5)
        m_idx = m_idx * lattice().size(5) + x5;
      if (ndim > 6)
        m_idx = m_idx * lattice().size(6) + x6;
      if (ndim > 7)
        m_idx = m_idx * lattice().size(7) + x7;
      if (ndim > 8)
        m_idx = m_idx * lattice().size(8) + x8;
      if (ndim > 9)
        m_idx = m_idx * lattice().size(9) + x9;

      m_idx = lattice().local(m_idx);

      if (m_idx == NOWHERE)
      {
        bool print = mdp.printing();
        mdp.enablePrinting();
        mdp << "Warning message from ME=" << ME << ":\n";
        mdp << "You assigned a site that is not here!\n";
        mdp.restorePrinting(print);
      }
    }

    bool operator==(const int *x)
    {
      mdp_int ndim = lattice().n_dimensions();
      bool is_it = true;
      for (mdp_int mu = 0; mu < ndim; mu++)
        if (x[mu] != lattice().coordinate(m_idx, mu))
          is_it = false;
      return is_it;
    }

    bool operator!=(const int *x)
    {
      return !(*this == x);
    }

    /** @brief checks the site coordinates vs the coordinates passed as args
     */
    int is_equal(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
                 int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
    {
      mdp_int ndim = lattice().n_dimensions();

      if ((ndim > 0) && (x0 != lattice().coordinate(m_idx, 0)))
        return false;
      if ((ndim > 1) && (x1 != lattice().coordinate(m_idx, 1)))
        return false;
      if ((ndim > 2) && (x2 != lattice().coordinate(m_idx, 2)))
        return false;
      if ((ndim > 3) && (x3 != lattice().coordinate(m_idx, 3)))
        return false;
      if ((ndim > 4) && (x4 != lattice().coordinate(m_idx, 4)))
        return false;
      if ((ndim > 5) && (x5 != lattice().coordinate(m_idx, 5)))
        return false;
      if ((ndim > 6) && (x6 != lattice().coordinate(m_idx, 6)))
        return false;
      if ((ndim > 7) && (x7 != lattice().coordinate(m_idx, 7)))
        return false;
      if ((ndim > 8) && (x8 != lattice().coordinate(m_idx, 8)))
        return false;
      if ((ndim > 9) && (x9 != lattice().coordinate(m_idx, 9)))
        return false;
      return true;
    }
  };

  /// converts a site into a binary number
  /// to be used only if the site is a vertex of an hypercube centered
  /// at the origin. this is used to make staggered mesons
  mdp_int site2binary(mdp_site x)
  {
    int a = 0;
    for (mdp_int mu = 0; mu < x.lattice().n_dimensions(); mu++)
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
  int on_which_process(mdp_lattice &a, int x0 = 0, int x1 = 0, int x2 = 0, int x3 = 0,
                       int x4 = 0, int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
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
    return a.where(x);
  }

#ifdef MDP_LATTICE

  /// Returns the local object mdp_prng at site x of the lattice
  mdp_prng &mdp_lattice::random(mdp_site x)
  {
    if (m_local_random_generator)
    {
      if (!x.is_in())
        error("request the random generator of a non local site");
      return m_random_obj[x.local_index() - start0(ME, 0)];
    }
    return mdp_random;
  }

#endif

  /// When compiled with TWISTED_BOUNDARY the mdp_site class keeps track of
  /// sites that moved around the boundary of the torus topology. this function
  /// returns false if this is one such site, true otherwise.
  int in_block(mdp_site x)
  {
#ifdef TWISTED_BOUNDARY
    for (mdp_int mu = 0; (mu < BLOCKSITE) && mu < x.lattice().n_dimensions(); mu++)
      if (x.block(mu) != 0)
        return false;
#endif
    return true;
  }

  std::ostream &operator<<(std::ostream &os, mdp_site &x)
  {
    for (mdp_int i = 0; i < x.lattice().n_dimensions(); i++)
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
