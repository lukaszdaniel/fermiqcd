/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class gauge_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_GAUGE_FIELD_
#define FERMIQCD_GAUGE_FIELD_

#include "mdp_complex_field.h"
#include "mdp_matrix_field.h"

namespace MDP
{
  /// @brief the chromo-electro-magnetic field for any SU(n)
  ///
  /// Example:
  /// @verbatim
  ///    int nc=3;
  ///    constexpr Box box = {10,8,8,8};
  ///    mdp_lattice lattice(box);
  ///    gauge_field U(lattice,nc);
  ///    mdp_site x(lattice);
  ///    U.load("myfield");
  ///    compute_em_field(U);
  ///    forallsites(x)
  ///      for(int mu=0; mu<U.ndim(); mu++)
  ///        for(int nu=mu+1; nu<U.ndim(); nu++)
  ///          std::cout << U.em(x,mu,nu) << std::endl;
  /// @endverbatim
  /// Note that U.em(x,mu,nu) is \f$ a^2 G_{\mu\nu} \f$ and
  /// it is a color matrix in SU(nc). \f$a\f$ is the lattice spacing.
  class em_field : public mdp_complex_field
  {
  private:
    mdp_int m_nc;

    int ordered_index(int mu, int nu) const
    {
      // ////////////////////////
      // this map mu, nu -> k  //
      //           0   1    0  //
      //           0   2    1  //
      //           0   3    2  //
      //           1   2    3  //
      //           1   3    4  //
      //           2   3    5  //
      // It works for any ndim //
      // ////////////////////////
      if (mu < nu)
        return nu + (mu * (2 * ndim() - mu - 3)) / 2 - 1;

      if (mu > nu)
      {
        mdp_int num_of_elements = ((ndim() * (ndim() - 1)) / 2);
        return mu + (nu * (2 * ndim() - nu - 3)) / 2 - 1 + num_of_elements;
      }
      // error("wrong call to ordered_index() with mu>=nu");
      return -1; // error in this case!
    }

  public:
    /** @brief Create a default em field
     */
    em_field() : mdp_complex_field(), m_nc(0)
    {
    }

    /** @brief Create an em field on \e nc_ colours
     *
     * @param a Lattice where the field needs to reside.
     * @param nc_ Number of colours.
     */
    em_field(mdp_lattice &a, int nc_) : mdp_complex_field(a, (nc_ * nc_ * ((a.ndim() * (a.ndim() - 1)) / 2))), m_nc(nc_)
    {
    }

    /** @brief Copy constructor
     *
     * @param em em field to be copied.
     */
    em_field(const em_field &em) : mdp_complex_field(em), m_nc(em.m_nc)
    {
    }

    void allocate_em_field(mdp_lattice &a, int nc_)
    {
      deallocate_field();
      m_nc = nc_;
      mdp_int num_of_elements = ((a.ndim() * (a.ndim() - 1)) / 2);
      allocate_field(a, num_of_elements * m_nc * m_nc);
    }

    mdp_int nc() const
    {
      return m_nc;
    }

    /** @brief returns the matrix in directions \e mu, \e nu stored at site x
     */
    mdp_matrix operator()(mdp_site x, int mu, int nu)
    {
#ifdef CHECK_ALL
      if (mu >= nu)
        error("em(x,mu,nu) for mu>=nu is not defined");
#endif
      int k = ordered_index(mu, nu);
      return mdp_matrix(address(x, k * m_nc * m_nc), m_nc, m_nc);
    }

    /** @brief returns the (i,j) component of the matrix in directions \e mu, \e nu stored at site x
     */
    mdp_complex &operator()(mdp_site x, int mu, int nu, int i, int j)
    {
#ifdef CHECK_ALL
      if (mu >= nu)
        error("em(x,mu,nu) for mu>=nu is not defined");
#endif
      int k = ordered_index(mu, nu);
      return *(address(x, (k * m_nc + i) * m_nc + j));
    }

    /** @brief returns the (i,j) const component of the matrix in directions \e mu, \e nu stored at site x
     */
    const mdp_complex &operator()(mdp_site x, int mu, int nu,
                                  int i, int j) const
    {
#ifdef CHECK_ALL
      if (mu >= nu)
        error("em(x,mu,nu) for mu>=nu is not defined");
#endif
      int k = ordered_index(mu, nu);
      return *(address(x, (k * m_nc + i) * m_nc + j));
    }
  };

#ifdef TWISTED_BOUNDARY
  void twist_boundary(mdp_matrix &M, mdp_site &x);
#endif

  /// @brief the gauge field for any SU(n)
  ///
  /// Example:
  /// @verbatim
  ///    int nc=3;
  ///    constexpr Box box = {10,8,8,8};
  ///    mdp_lattice lattice(box);
  ///    gauge_field U(lattice,nc);
  ///    mdp_site x(lattice);
  ///    // set_cold(U);
  ///    forallsites(x)
  ///       for(int mu=0; mu<U.ndim(); mu++)
  ///          U(x,mu)=1;
  ///    U.update(); // synchronization
  ///    U.save("myfield");
  ///    U.load("myfield");
  /// @endverbatim
  /// Note that U(x,mu) is \f$ \exp{iaA_{\mu}} \f$ and
  /// it is a color matrix in SU(nc). \f$a\f$ is the lattice spacing.
  class gauge_field : public mdp_complex_field
  {
  private:
    mdp_int m_nc;

  public:
    em_field em;
    mdp_nmatrix_field long_links;
    mdp_int_scalar_field i_jump;
    mdp_matrix_field swirls;

    /** @brief Create a default gauge field
     */
    gauge_field() : mdp_complex_field(), m_nc(0)
    {
    }

    /** @brief Create a gauge field on \e nc_ colours
     *
     * @param a Lattice where the field needs to reside.
     * @param nc_ Number of colours.
     */
    gauge_field(mdp_lattice &a, int nc_) : mdp_complex_field(a, (nc_ * nc_ * a.ndim())), m_nc(nc_)
    {
    }

    /** @brief Copy constructor
     *
     * @param U Gauge field to be copied.
     */
    gauge_field(const gauge_field &U) : mdp_complex_field(U), m_nc(U.m_nc)
    {
    }

    gauge_field &operator=(const gauge_field &U)
    {
      if (this == &U)
        return *this;

      // base assignment
      mdp_complex_field::operator=(U);

      // derived fields assignment
      m_nc = U.m_nc;

      return *this;
    }

    void allocate_gauge_field(mdp_lattice &a, int nc_)
    {
      deallocate_field();
      m_nc = nc_;
      allocate_field(a, a.ndim() * m_nc * m_nc);
    }

    mdp_int nc() const
    {
      return m_nc;
    }

    /** @brief returns the matrix in direction \e mu stored at site x
     */
    mdp_matrix operator()(mdp_site x, int mu)
    {
#ifndef TWISTED_BOUNDARY
      return mdp_matrix(address(x, mu * m_nc * m_nc), m_nc, m_nc);
#else
      mdp_matrix tmp(address(x, mu * m_nc * m_nc), m_nc, m_nc);
      if (!in_block(x))
      {
        mdp_matrix a;
        a = tmp;
#ifdef ADVANCED_TWIST
        twist_eat_links(a, x, mu);
#else
        twist_boundary(a, x);
#endif
        return (a);
      }
      return tmp;
#endif
    }

    /** @brief returns the const matrix in direction \e mu stored at site x
     */
    const mdp_matrix operator()(mdp_site x, int mu) const
    {
#ifndef TWISTED_BOUNDARY
      return mdp_matrix(address(x, mu * m_nc * m_nc), m_nc, m_nc);
#else
      mdp_matrix tmp(address(x, mu * m_nc * m_nc), m_nc, m_nc);
      if (!in_block(x))
      {
        mdp_matrix a;
        a = tmp;
#ifdef ADVANCED_TWIST
        twist_eat_links(a, x, mu);
#else
        twist_boundary(a, x);
#endif
        return (a);
      }
      return tmp;
#endif
    }

    /** @brief returns the (i,j) component of the matrix in direction \e mu stored at site x
     */
    mdp_complex &operator()(mdp_site x, int mu, int i, int j)
    {
#ifdef TWISTED_BOUNDARY
      if (!in_block(x))
        error("call to &operator=(mdp_site x) per x out of block");
#endif
      return *(address(x, (mu * m_nc + i) * m_nc + j));
    }

    /** @brief returns the (i,j) const component of the matrix in direction \e mu stored at site x
     */
    const mdp_complex &operator()(mdp_site x, int mu, int i, int j) const
    {
#ifdef TWISTED_BOUNDARY
      if (!in_block(x))
        error("call to &operator=(mdp_site x) per x out of block");
#endif
      return *(address(x, (mu * m_nc + i) * m_nc + j));
    }

    /** @brief returns the matrix in direction \e mu stored at site x
     *
     * @note if \e sign is negative returned matrix is a hermitian matrix
     * if direction \e -mu
     */
    mdp_matrix operator()(mdp_site x, int sign, int mu)
    {
      if (sign == +1)
        return (*this)(x, mu);
      if (sign == -1)
        return hermitian((*this)(x - mu, mu));

      return mdp_matrix();
    }

    /** @brief returns the const matrix in direction \e mu stored at site x
     *
     * @note if \e sign is negative returned matrix is a hermitian matrix
     * if direction \e -mu
     */
    const mdp_matrix operator()(mdp_site x, int sign, int mu) const
    {
      if (sign == +1)
        return (*this)(x, mu);
      if (sign == -1)
        return hermitian((*this)(x - mu, mu));

      return mdp_matrix();
    }

#ifndef TWISTED_BOUNDARY
    /** @brief returns the (i,j) const component of the matrix in direction \e mu stored at site x
     *
     * @note if \e sign is negative returned element is conjugated
     */
    const mdp_complex operator()(mdp_site x, int sign, int mu,
                                 int i, int j) const
    {
      if (sign == +1)
        return *(address(x, (mu * m_nc + i) * m_nc + j));
      if (sign == -1)
        return conj(*(address(x - mu, (mu * m_nc + j) * m_nc + i)));

      error("call to U(x,0,mu,i,j)");

      return mdp_complex(0, 0);
    }
#endif
  };

  // /////////////////////////////////////////////////////
  // /////////////////////////////////////////////////////
  // stuff for twisted boundary conditions
  // ignore if you do not care
  // /////////////////////////////////////////////////////
  // /////////////////////////////////////////////////////

#ifdef TWISTED_BOUNDARY
  void twist_boundary(mdp_matrix &M, mdp_site &x)
  {
    static int block;
    static mdp_complex z = exp(mdp_complex(0, 2.0 * Pi / 3.0));
    static mdp_complex a, b, c, d, e, f, g, h, i;

    for (mdp_int mu = 1; mu < x.lattice().ndim(); mu++)
    {
      block = x.block(mu);
      if (block != 0)
      {
        a = M(0, 0);
        b = M(0, 1);
        c = M(0, 2);
        d = M(1, 0);
        e = M(1, 1);
        f = M(1, 2);
        g = M(2, 0);
        h = M(2, 1);
        i = M(2, 2);
        if (mu == 1 && block == 1)
        {
          M(0, 0) = e;
          M(0, 1) = f;
          M(0, 2) = d;
          M(1, 0) = h;
          M(1, 1) = i;
          M(1, 2) = g;
          M(2, 0) = b;
          M(2, 1) = c;
          M(2, 2) = a;
        }
        if (mu == 1 && block == -1)
        {
          M(0, 0) = i;
          M(0, 1) = g;
          M(0, 2) = h;
          M(1, 0) = c;
          M(1, 1) = a;
          M(1, 2) = b;
          M(2, 0) = f;
          M(2, 1) = d;
          M(2, 2) = e;
        }
        if (mu == 2 && block == 1)
        {
          M(0, 0) = a;
          M(0, 1) = b / z;
          M(0, 2) = c / z / z;
          M(1, 0) = d * z;
          M(1, 1) = e;
          M(1, 2) = f / z;
          M(2, 0) = g * z * z;
          M(2, 1) = h * z;
          M(2, 2) = i;
        }
        if (mu == 2 && block == -1)
        {
          M(0, 0) = a;
          M(0, 1) = b * z;
          M(0, 2) = c * z * z;
          M(1, 0) = d / z;
          M(1, 1) = e;
          M(1, 2) = f * z;
          M(2, 0) = g / z / z;
          M(2, 1) = h / z;
          M(2, 2) = i;
        }
        if (mu == 3 && block == 1)
        {
          M(0, 0) = i;
          M(0, 1) = g / z;
          M(0, 2) = h / z / z;
          M(1, 0) = c * z;
          M(1, 1) = a;
          M(1, 2) = b / z;
          M(2, 0) = f * z * z;
          M(2, 1) = d * z;
          M(2, 2) = e;
        }
        if (mu == 3 && block == -1)
        {
          M(0, 0) = e;
          M(0, 1) = f * z;
          M(0, 2) = d * z * z;
          M(1, 0) = h / z;
          M(1, 1) = i;
          M(1, 2) = g * z;
          M(2, 0) = b / z / z;
          M(2, 1) = c / z;
          M(2, 2) = a;
        }
        if (block * block > 1)
          error("two blocks out\n");
      }
    }
  }

  void define_twist_matrices()
  {
    begin_function("gauge_field__define_twist_matrices");
    mdp_complex z = exp(mdp_complex(0, 2.0 * Pi / 3.0));
    OmegaTwist[0] = mdp_identity(3);
    OmegaTwist[1] = mdp_zero(3);
    OmegaTwist[1](0, 1) = 1;
    OmegaTwist[1](1, 2) = 1;
    OmegaTwist[1](2, 0) = 1;
    OmegaTwist[2] = mdp_zero(3);
    OmegaTwist[2](0, 0) = 1.0 / z;
    OmegaTwist[2](1, 1) = 1;
    OmegaTwist[2](2, 2) = z;
    OmegaTwist[3] = mdp_zero(3);
    OmegaTwist[3](0, 2) = 1.0 / z;
    OmegaTwist[3](1, 0) = 1;
    OmegaTwist[3](2, 1) = z;
    end_function("gauge_field__define_twist_matrices");
  }

  void twist_eat_fields(mdp_matrix &M, mdp_site &x,
                        gauge_field &omega)
  {
    begin_function("gauge_field__twist_eat_matrices");
    static int block;
    // static mdp_complex z = exp(mdp_complex(0, 2.0 * Pi / 3.0));
    static mdp_complex a, b, c, d, e, f, g, h, i;

    for (mdp_int mu = 1; mu < x.lattice().ndim(); mu++)
    {
      block = x.block(mu);
      if (block != 0)
      {
        if (mu == 1 && block == 1)
          M = omega(x, mu) * M * hermitian(omega(x, mu));
        if (mu == 1 && block == -1)
          M = hermitian(omega(x, mu)) * M * omega(x, mu);
        if (mu == 2 && block == 1)
          M = omega(x, mu) * M * hermitian(omega(x, mu));
        if (mu == 2 && block == -1)
          M = hermitian(omega(x, mu)) * M * omega(x, mu);
        if (mu == 3 && block == 1)
          M = omega(x, mu) * M * hermitian(omega(x, mu));
        if (mu == 3 && block == -1)
          M = hermitian(omega(x, mu)) * M * omega(x, mu);
        if (block * block > 1)
          error("two blocks out\n");
      }
    }
    end_function("gauge_field__twist_eat_matrices");
  }
#endif
} // namespace MDP

#endif /* FERMIQCD_GAUGE_FIELD_ */
