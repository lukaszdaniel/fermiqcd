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
    mdp_real &real() { return m_re; }
    mdp_real &imag() { return m_im; }

    const mdp_real &real() const { return m_re; }
    const mdp_real &imag() const { return m_im; }

    /** @brief Create a complex number
     *
     * @param a Real part of the complex number
     * @param b Imaginary part of the complex number
     */
    mdp_complex(const mdp_real a = 0.0, const mdp_real b = 0.0) : m_re(a), m_im(b)
    {
    }

    /** @brief Copy constructor
     */
    mdp_complex(const mdp_complex &c) : m_re(c.m_re), m_im(c.m_im)
    {
    }

    bool operator==(const mdp_complex &c) const
    {
      return ((m_re == c.m_re) && (m_im == c.m_im));
    }

    bool operator!=(const mdp_complex &c) const
    {
      return !(*this == c);
    }

    mdp_complex operator+() const
    {
      return (*this);
    }

    mdp_complex operator-() const
    {
      return mdp_complex(-m_re, -m_im);
    }

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
      (*this) = (*this) * c;
    }

    void operator/=(const mdp_complex &c)
    {
      (*this) = (*this) / c;
    }
  };

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
    mdp_real r = std::pow((double)abs(c), (double)z), p = (mdp_real)z * phase(c);
    return mdp_complex(r * std::cos(p), r * std::sin(p));
  }

  mdp_complex sqrt(const mdp_complex &c)
  {
    mdp_real r = std::sqrt(abs(c)), p = 0.5 * arg(c);
    return mdp_complex(r * std::cos(p), r * std::sin(p));
  }

  mdp_complex exp(const mdp_complex &c)
  {
    mdp_real exp_c_re = std::exp(c.real());
    return mdp_complex(exp_c_re * std::cos(c.imag()), exp_c_re * std::sin(c.imag()));
  }

  mdp_complex sin(const mdp_complex &c)
  {
    return mdp_complex(cosh(c.imag()) * std::sin(c.real()), sinh(c.imag()) * std::cos(c.real()));
  }

  mdp_complex cos(const mdp_complex &c)
  {
    return mdp_complex(cosh(c.imag()) * std::cos(c.real()), -sinh(c.imag()) * std::sin(c.real()));
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
#endif

  mdp_complex operator+(const mdp_complex &a, const int c)
  {
    return mdp_complex(a.real() + c, a.imag());
  }

  mdp_complex operator-(const mdp_complex &a, const int c)
  {
    return mdp_complex(a.real() - c, a.imag());
  }

  mdp_complex operator*(const mdp_complex &a, const int c)
  {
    return mdp_complex(a.real() * c, a.imag() * c);
  }

  mdp_complex operator/(const mdp_complex &a, const int c)
  {
    return mdp_complex(a.real() / c, a.imag() / c);
  }

  mdp_complex operator+(const int a, const mdp_complex &c)
  {
    return mdp_complex(c.real() + a, c.imag());
  }

  mdp_complex operator-(const int a, const mdp_complex &c)
  {
    return mdp_complex(a - c.real(), -c.imag());
  }

  mdp_complex operator*(const int a, const mdp_complex &c)
  {
    return mdp_complex(c.real() * a, c.imag() * a);
  }

  mdp_complex operator/(const int a, const mdp_complex &c)
  {
    mdp_real den = (c.real() * c.real() + c.imag() * c.imag()) / a;
    return mdp_complex((c.real()) / den, (-c.imag()) / den);
  }

  mdp_complex operator+(const mdp_complex &a, const float c)
  {
    return mdp_complex(a.real() + c, a.imag());
  }

  mdp_complex operator-(const mdp_complex &a, const float c)
  {
    return mdp_complex(a.real() - c, a.imag());
  }

  mdp_complex operator*(const mdp_complex &a, const float c)
  {
    return mdp_complex(a.real() * c, a.imag() * c);
  }

  mdp_complex operator/(const mdp_complex &a, const float c)
  {
    return mdp_complex(a.real() / c, a.imag() / c);
  }

  mdp_complex operator+(const float a, const mdp_complex &c)
  {
    return mdp_complex(c.real() + a, c.imag());
  }

  mdp_complex operator-(const float a, const mdp_complex &c)
  {
    return mdp_complex(a - c.real(), -c.imag());
  }

  mdp_complex operator*(const float a, const mdp_complex &c)
  {
    return mdp_complex(c.real() * a, c.imag() * a);
  }

  mdp_complex operator/(const float a, const mdp_complex &c)
  {
    mdp_real den = (c.real() * c.real() + c.imag() * c.imag()) / a;
    return mdp_complex((c.real()) / den, (-c.imag()) / den);
  }

  mdp_complex operator+(const mdp_complex &a, const double c)
  {
    return mdp_complex(a.real() + c, a.imag());
  }

  mdp_complex operator-(const mdp_complex &a, const double c)
  {
    return mdp_complex(a.real() - c, a.imag());
  }

  mdp_complex operator*(const mdp_complex &a, const double c)
  {
    return mdp_complex(a.real() * c, a.imag() * c);
  }

  mdp_complex operator/(const mdp_complex &a, const double c)
  {
    return mdp_complex(a.real() / c, a.imag() / c);
  }

  mdp_complex operator+(const double a, const mdp_complex &c)
  {
    return mdp_complex(c.real() + a, c.imag());
  }

  mdp_complex operator-(const double a, const mdp_complex &c)
  {
    return mdp_complex(a - c.real(), -c.imag());
  }

  mdp_complex operator*(const double a, const mdp_complex &c)
  {
    return mdp_complex(c.real() * a, c.imag() * a);
  }

  mdp_complex operator/(const double a, const mdp_complex &c)
  {
    mdp_real den = (c.real() * c.real() + c.imag() * c.imag()) / a;
    return mdp_complex((c.real()) / den, (-c.imag()) / den);
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
      os << a.real() << "+" << std::abs(a.imag()) << "I";
    return os;
  }
} // namespace MDP

#endif /* MDP_COMPLEX_ */
