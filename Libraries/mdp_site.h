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

#include <span>
#include "mdp_vector.h"
#include "mdp_lattice.h"

namespace MDP
{
  /// @brief site object to loop on a lattice
  ///
  /// Example:
  /// @verbatim
  ///   constexpr Box box = {10,10,10};
  ///   mdp_lattice lattice(box);
  ///   mdp_site x(lattice);
  ///   forallsites(x) std::cout << x << std::endl;
  ///   if (on_which_process(lattice,1,1,1)==ME) {
  ///      x.set(1,1,1);
  ///      std::cout << lattice.random(x).plain() << std::endl;
  ///   }
  /// @endverbatim
  class mdp_site
  {
  private:
    const mdp_lattice *m_lattice; /// this points to the lattice for this field
    mdp_int m_idx;                /// value of the mdp_site in local coordinate

#ifdef BLOCKSITE
    int m_block[BLOCKSITE];
#endif

    /** @brief Coordinte setter for mdp_site
     *
     * @tparam Accessor type of the coordinate accessor, must be callable with signature: mdp_int(size_t i)
     *
     * This function sets the site to a the location specified by the coordinates provided by the accessor.
     * The accessor is called with indices from 0 to (ndim-1) to retrieve the coordinates in each dimension.
     */
    template <typename Accessor>
    void set_from_accessor(Accessor coord)
    {
      const auto &lat = lattice();
      const mdp_int ndim = lat.n_dimensions();

      m_idx = coord(0);

      for (mdp_int i = 1; i < ndim; ++i)
        m_idx = m_idx * lat.size(i) + coord(i);

      m_idx = lat.local(m_idx);

      if (m_idx == NOWHERE)
      {
        bool print = mdp.printingEnabled();
        mdp.enablePrinting();
        mdp << "Warning message from ME=" << ME << ":\n";
        mdp << "You assigned a site that is not here!\n";
        mdp.setPrinting(print);
      }
    }

    /** @brief checks if the site is equal to a site defined by coordinates
     *
     * @tparam Accessor type of the coordinate accessor, must be callable with signature: mdp_int(size_t i)
     */
    template <typename Accessor>
    bool is_equal_from_accessor(Accessor coord) const
    {
      const auto &lat = lattice();
      const mdp_int ndim = lat.n_dimensions();

      for (mdp_int i = 0; i < ndim; ++i)
        if (coord(i) != lat.coordinate(m_idx, i))
          return false;

      return true;
    }

  public:
    /** @brief declares object of class mdp_site living on the
     * lattice passed by reference
     */
    mdp_site(const mdp_lattice &a)
    {
      m_lattice = &a;
      m_idx = (*m_lattice).start0(ME, 0);
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = 0;
#endif
    }

    mdp_site(mdp_int i, const mdp_lattice *ptr2)
    {
      m_idx = i;
      m_lattice = ptr2;
#ifdef BLOCKSITE
      for (int k = 0; k < BLOCKSITE; k++)
        m_block[k] = 0;
#endif
    }

#ifdef BLOCKSITE
    mdp_site(mdp_int i, const mdp_lattice *ptr2, int b[], int sign = 0, int mu = 0)
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
    const mdp_lattice &lattice() const
    {
      return *m_lattice;
    }

#ifdef BLOCKSITE
    int block(mdp_int mu) const
    {
      if (mu >= BLOCKSITE)
        error("BLOCKSITE < lattice dimension");

      return m_block[mu];
    }
#endif

    mdp_site &operator=(mdp_int i)
    {
      m_idx = lattice().start0(ME, 0) + i;
      return *this;
    }

    mdp_site &operator=(const mdp_site &x)
    {
      if (this == &x)
        return *this;

      if (m_lattice == x.m_lattice)
        m_idx = x.m_idx;
      else
        set_global(x.global_index());

#ifdef BLOCKSITE
      std::copy(std::begin(x.m_block), std::end(x.m_block), m_block);
#endif

      return *this;
    }

    bool operator==(const mdp_site &x) const
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
    bool is_in() const
    {
      return ((m_idx >= lattice().start0(ME, 0)) && (m_idx < lattice().stop0(ME, 1)));
    }

    /** @brief checks if the site is inside the portion of the lattice stored by
     * the current process or if the site is in a local copy of a remote
     * site
     */
    int is_here() const
    {
      return ((m_idx >= 0) && (m_idx < lattice().enclosing_volume()));
    }

    /** @brief returns the parity EVEN or ODD of the site
     */
    mdp_int parity() const
    {
      return lattice().site_parity(m_idx);
    }

    /** @brief true if the site is stored locally as a copy of a site local
     * in another process
     */
    bool is_in_boundary() const
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
    mdp_int global_index() const
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
     * in direction mu=(0...ndim-1)
     */
    mdp_site hop(int n, mdp_int mu) const
    {
      mdp_site y = *this;

      if (n > 0)
        while (n--)
          y = y + mu;
      else
        while (n++)
          y = y - mu;

      return y;
    }

    /** @brief sets the site to the coordinates stored in vector v
     */
    mdp_site &operator=(const mdp_vector &v)
    {
      set(v);
      return *this;
    }

    /** @brief returns a site similar to the present but
     * each coordinates mu of the site shifted according to v[mu]
     */
    mdp_site operator+(const mdp_vector &v)
    {
      mdp_site y = *this;
      for (mdp_uint mu = 0; mu < lattice().n_dimensions(); ++mu)
      {
        for (mdp_int step = 0; step < std::abs(v[mu]); ++step)
          y = (v[mu] > 0) ? (y + mu) : (y - mu);
      }
      return y;
    }

    /** @brief returns a site similar to the present but
     * each coordinates mu of the site shifted according to -v[mu]
     */
    mdp_site operator-(const mdp_vector &v)
    {
      mdp_site y = *this;
      for (mdp_uint mu = 0; mu < lattice().n_dimensions(); ++mu)
      {
        for (mdp_int step = 0; step < std::abs(v[mu]); ++step)
          y = (v[mu] > 0) ? (y - mu) : (y + mu);
      }
      return y;
    }

    /// returns mu coordinate of the site
    mdp_int operator()(mdp_int mu) const
    {
      return lattice().coordinate(m_idx, mu);
    }

    mdp_site &operator=(std::span<const int> x)
    {
      set_from_accessor([&](mdp_int i)
                        { return x[i]; });
      return *this;
    }

    void set(const mdp_vector &v)
    {
      set_from_accessor([&](mdp_int i)
                        { return v[i]; });
    }

    /// sets the site to a the location specified by the coordinates
    /// and assumes the site is local (or at least a copy)
    /// @see on_which_process()
    template <typename... Args>
    void set(Args... args)
    {
      static_assert(sizeof...(Args) <= mdp_vector::dim);
      mdp_vector v{args...};
      set(v);
    }

    bool operator==(const int x[])
    {
      mdp_int ndim = lattice().n_dimensions();
      bool is_it = true;
      for (mdp_int mu = 0; mu < ndim; mu++)
        if (x[mu] != lattice().coordinate(m_idx, mu))
          is_it = false;
      return is_it;
    }

    bool operator!=(const int x[])
    {
      return !(*this == x);
    }

    /** @brief checks the site coordinates vs the coordinates passed as args
     */
    bool is_equal(const mdp_vector &v) const
    {
      return is_equal_from_accessor([&](mdp_int i)
                                    { return v[i]; });
    }

    bool is_equal(std::span<const int> x) const
    {
      return is_equal_from_accessor([&](mdp_int i)
                                    { return x[i]; });
    }

    template <typename... Args>
    bool is_equal(Args... args) const
    {
      static_assert(sizeof...(Args) > 0);
      static_assert(sizeof...(Args) <= VECTOR_MAX_DIM);

      std::array<int, sizeof...(Args)> coords{args...};
      return is_equal(std::span<const int>(coords));
    }
  };

  /// converts a site into a binary number
  /// to be used only if the site is a vertex of an hypercube centered
  /// at the origin. this is used to make staggered mesons
  mdp_int site2binary(mdp_site x)
  {
    int a = 0;
    for (mdp_uint mu = 0; mu < x.lattice().n_dimensions(); mu++)
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
  int on_which_process(const mdp_lattice &a, std::span<const int> x)
  {
    return a.where(x.data());
  }

  template <typename... Args>
  int on_which_process(const mdp_lattice &a, Args... args)
  {
    static_assert(sizeof...(Args) > 0);
    static_assert(sizeof...(Args) <= VECTOR_MAX_DIM);

    std::array<int, sizeof...(Args)> coords{args...};
    return on_which_process(a, std::span<const int>(coords));
  }

#ifdef MDP_LATTICE

  /// Returns the local object mdp_random at site x of the lattice
  mdp_random &mdp_lattice::random(mdp_site x) const
  {
    if (m_local_random_generator)
    {
      if (!x.is_in())
        error("request the random generator of a non local site");
      return m_random_obj[x.local_index() - start0(ME, 0)];
    }
    return mdp_global_random;
  }

#endif

  /// When compiled with TWISTED_BOUNDARY the mdp_site class keeps track of
  /// sites that moved around the boundary of the torus topology. this function
  /// returns false if this is one such site, true otherwise.
#ifdef TWISTED_BOUNDARY
  int in_block(mdp_site x)
  {
    for (mdp_uint mu = 0; (mu < BLOCKSITE) && mu < x.lattice().n_dimensions(); mu++)
      if (x.block(mu) != 0)
        return false;
    return true;
  }
#endif

  std::ostream &operator<<(std::ostream &os, mdp_site &x)
  {
    for (mdp_uint i = 0; i < x.lattice().n_dimensions(); i++)
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
