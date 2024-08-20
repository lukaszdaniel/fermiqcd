/////////////////////////////////////////////////////////////////
/// @file mdp_complex.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains delcaration of class mdp_complex for complex numbers
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_COMPLEX_
#define MDP_COMPLEX_

// #define DO_NOT_USE_MDP_COMPLEX

#include <ostream>
#include <cmath>

#ifdef DO_NOT_USE_MDP_COMPLEX
#include <complex>
#endif
#include "mdp_global_vars.h"

namespace MDP
{
#ifdef DO_NOT_USE_MDP_COMPLEX
  typedef std::complex<mdp_real> mdp_complex;
#else

  /// @brief portable complex numbers
  ///
  /// Example:
  /// @verbatim
  ///    mdp_complex x=3+5*I;
  ///    cout << x.real() << "," << x.imag() << endl;
  ///    cout << sin(x) << endl;
  /// @endverbatim
  class mdp_complex
  {
  private:
    mdp_real m_re;
    mdp_real m_im;

  public:
    // Accessors for real and imaginary parts
    mdp_real &real() { return m_re; }
    mdp_real &imag() { return m_im; }
    const mdp_real &real() const { return m_re; }
    const mdp_real &imag() const { return m_im; }

    /** @brief Create a complex number
     *
     * @param a Real part of the complex number
     * @param b Imaginary part of the complex number
     */
    mdp_complex(mdp_real a = 0.0, mdp_real b = 0.0) : m_re(a), m_im(b) {}

    /** @brief Copy constructor
     */
    mdp_complex(const mdp_complex &c) : m_re(c.m_re), m_im(c.m_im) {}

    // Comparison operators
    bool operator==(const mdp_complex &c) const
    {
      return (m_re == c.m_re) && (m_im == c.m_im);
    }

    bool operator!=(const mdp_complex &c) const
    {
      return !(*this == c);
    }

    // Unary operators
    mdp_complex operator+() const { return *this; }
    mdp_complex operator-() const { return mdp_complex(-m_re, -m_im); }

    // Arithmetic operators
    mdp_complex operator+(const mdp_complex &a) const
    {
      return mdp_complex(m_re + a.m_re, m_im + a.m_im);
    }

    mdp_complex operator-(const mdp_complex &a) const
    {
      return mdp_complex(m_re - a.m_re, m_im - a.m_im);
    }

    mdp_complex operator*(const mdp_complex &c) const
    {
      mdp_real re = m_re * c.m_re - m_im * c.m_im;
      mdp_real im = m_re * c.m_im + m_im * c.m_re;
      return mdp_complex(re, im);
    }

    mdp_complex operator/(const mdp_complex &c) const
    {
      mdp_real norm2 = c.m_re * c.m_re + c.m_im * c.m_im;
      mdp_real re = (m_re * c.m_re + m_im * c.m_im) / norm2;
      mdp_real im = (m_im * c.m_re - m_re * c.m_im) / norm2;
      return mdp_complex(re, im);
    }

    // Compound assignment operators
    void operator+=(const mdp_complex &c)
    {
      m_re += c.m_re;
      m_im += c.m_im;
    }

    void operator-=(const mdp_complex &c)
    {
      m_re -= c.m_re;
      m_im -= c.m_im;
    }

    void operator*=(const mdp_complex &c)
    {
      *this = *this * c;
    }

    void operator/=(const mdp_complex &c)
    {
      *this = *this / c;
    }
  };

  // Helper functions for mdp_complex
  mdp_real real(const mdp_complex &c)
  {
    return c.real();
  }

  mdp_real imag(const mdp_complex &c)
  {
    return c.imag();
  }

  mdp_real abs(const mdp_complex &c)
  {
    return std::sqrt(c.real() * c.real() + c.imag() * c.imag());
  }

  mdp_real phase(const mdp_complex &c)
  {
    return atan2(c.imag(), c.real());
  }

  mdp_real arg(const mdp_complex &c)
  {
    return phase(c);
  }

  mdp_complex pow(const mdp_complex &c, mdp_real z)
  {
    mdp_real r = std::pow(abs(c), z);
    mdp_real p = z * phase(c);
    return mdp_complex(r * std::cos(p), r * std::sin(p));
  }

  mdp_complex sqrt(const mdp_complex &c)
  {
    mdp_real r = std::sqrt(abs(c));
    mdp_real p = 0.5 * arg(c);
    return mdp_complex(r * std::cos(p), r * std::sin(p));
  }

  mdp_complex exp(const mdp_complex &c)
  {
    mdp_real exp_re = std::exp(c.real());
    return mdp_complex(exp_re * std::cos(c.imag()), exp_re * std::sin(c.imag()));
  }

  mdp_complex sin(const mdp_complex &c)
  {
    return mdp_complex(std::cosh(c.imag()) * std::sin(c.real()), std::sinh(c.imag()) * std::cos(c.real()));
  }

  mdp_complex cos(const mdp_complex &c)
  {
    return mdp_complex(std::cosh(c.imag()) * std::cos(c.real()), -std::sinh(c.imag()) * std::sin(c.real()));
  }

  mdp_complex times_i(const mdp_complex &c)
  {
    return mdp_complex(-c.imag(), c.real());
  }

  mdp_complex times_minus_i(const mdp_complex &c)
  {
    return mdp_complex(c.imag(), -c.real());
  }

  mdp_complex conj(const mdp_complex &a)
  {
    return mdp_complex(a.real(), -a.imag());
  }
#endif // DO_NOT_USE_MDP_COMPLEX

  // Operator overloads for mdp_complex and scalar (double) arithmetic
  mdp_complex operator+(const mdp_complex &a, double c)
  {
    return mdp_complex(a.real() + c, a.imag());
  }

  mdp_complex operator+(double a, const mdp_complex &c)
  {
    return mdp_complex(a + c.real(), c.imag());
  }

  mdp_complex operator-(const mdp_complex &a, double c)
  {
    return mdp_complex(a.real() - c, a.imag());
  }

  mdp_complex operator-(double a, const mdp_complex &c)
  {
    return mdp_complex(a - c.real(), -c.imag());
  }

  mdp_complex operator*(const mdp_complex &a, double c)
  {
    return mdp_complex(a.real() * c, a.imag() * c);
  }

  mdp_complex operator*(double a, const mdp_complex &c)
  {
    return mdp_complex(a * c.real(), a * c.imag());
  }

  mdp_complex operator/(const mdp_complex &a, double c)
  {
    return mdp_complex(a.real() / c, a.imag() / c);
  }

  mdp_complex operator/(double a, const mdp_complex &c)
  {
    double norm2 = c.real() * c.real() + c.imag() * c.imag();
    double re = a * c.real() / norm2;
    double im = -a * c.imag() / norm2;
    return mdp_complex(re, im);
  }

  const mdp_complex I = mdp_complex(0, 1);

  mdp_real abs2(const mdp_complex &a)
  {
    return a.real() * a.real() + a.imag() * a.imag();
  }

  std::ostream &operator<<(std::ostream &os, const mdp_complex &a)
  {
    if (a.imag() < 0)
      os << a.real() << "-" << -a.imag() << "I";
    else
      os << a.real() << "+" << a.imag() << "I";
    return os;
  }
} // namespace MDP

#endif /* MDP_COMPLEX_ */
