/////////////////////////////////////////////////////////////////
/// @file mdp_measure.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_measure
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_MEASURE_
#define MDP_MEASURE_

#include <cmath>
#include <iostream>
#include "mdp_compatibility_macros.h"
#include "mdp_random.h"

namespace MDP
{
  /// @brief implements error propagation
  ///
  /// Example:
  /// @verbatim
  ///    mdp_measure m;
  ///    // store 10 measurements
  ///    for(mdp_uint i=0; i<10; i++)
  ///       m << 3.0+mdp_global_random.gaussian(2.0);
  ///    m=sin(exp(m)+m);
  ///    std::cout << m.getmean() << "+/-" << m.getmerr() << std::endl;
  /// @endverbatim
  /// Assumes gaussian error propagation
  class mdp_measure
  {
  private:
    mdp_real m_mean;
    mdp_real m_error;
    mdp_uint m_num;

  public:
    mdp_measure(mdp_real mean_ = 0, mdp_real error_ = 0, mdp_uint num_ = 0) : m_mean(mean_), m_error(error_), m_num(num_)
    {
    }

    mdp_uint getnum() const
    {
      return m_num;
    }

    mdp_real getmean() const
    {
      return m_mean;
    }

    void setmean(mdp_real mean)
    {
      m_mean = mean;
    }

    mdp_real getmerr() const
    {
      return m_error;
    }

    void setmerror(mdp_real err)
    {
      m_error = err;
    }

    void reset()
    {
      m_num = 0;
      m_mean = 0;
      m_error = 0;
    }

    void set(mdp_real x, mdp_real dx, mdp_uint i = 1)
    {
      m_num = i;
      m_mean = x;
      m_error = dx;
    }

    void operator<<(mdp_real x)
    {
      mdp_real err2 = m_num * (std::pow(m_error, 2.0) + m_mean * m_mean) + std::pow(x, 2.0);
      ++m_num;
      m_mean = (m_mean * (m_num - 1) + x) / m_num;
      m_error = std::sqrt(err2 / m_num - m_mean * m_mean);
    }

    void operator>>(mdp_real &x)
    {
      x = m_mean + m_error * Random.gaussian();
    }

    mdp_measure operator+()
    {
      return mdp_measure(m_mean, m_error);
    }

    mdp_measure operator-()
    {
      return mdp_measure(-m_mean, m_error);
    }

    mdp_measure operator+(mdp_measure b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean + b.m_mean;
      tmp.m_error = m_error + b.m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator-(mdp_measure b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean - b.m_mean;
      tmp.m_error = m_error + b.m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator*(mdp_measure b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean * b.m_mean;
      tmp.m_error = std::fabs(m_mean) * b.m_error + m_error * std::fabs(b.m_mean);
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator/(mdp_measure b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean / b.m_mean;
      tmp.m_error = m_error / std::fabs(b.m_mean) +
                    std::fabs(m_mean) / std::pow(b.m_mean, 2.0) * b.m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator+(mdp_real b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean + b;
      tmp.m_error = m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator-(mdp_real b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean - b;
      tmp.m_error = m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator*(mdp_real b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean * b;
      tmp.m_error = m_error * std::fabs(b);
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator/(mdp_real b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean / b;
      tmp.m_error = m_error / std::fabs(b);
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure exp()
    {
      mdp_measure tmp;
      tmp.m_mean = std::exp(m_mean);
      tmp.m_error = std::exp(m_mean) * m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure log()
    {
      mdp_measure tmp;
      tmp.m_mean = std::log(m_mean);
      tmp.m_error = m_error / fabs(m_mean);
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure sin()
    {
      mdp_measure tmp;
      tmp.m_mean = std::sin(m_mean);
      tmp.m_error = std::fabs(std::cos(m_mean)) * m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure cos()
    {
      mdp_measure tmp;
      tmp.m_mean = std::cos(m_mean);
      tmp.m_error = std::fabs(std::sin(m_mean)) * m_error;
      tmp.m_num = 1;
      return tmp;
    }
  };

  mdp_measure operator+(mdp_real a, mdp_measure b)
  {
    return b + a;
  }

  mdp_measure operator-(mdp_real a, mdp_measure b)
  {
    return -(b - a);
  }

  mdp_measure operator*(mdp_real a, mdp_measure b)
  {
    return b * a;
  }

  mdp_measure operator/(mdp_real a, mdp_measure b)
  {
    mdp_real mean = a / b.getmean();
    mdp_real error = std::fabs(a) / std::pow(b.getmean(), 2.0) * b.getmerr();
    return mdp_measure(mean, error);
  }

  mdp_measure exp(mdp_measure a)
  {
    return a.exp();
  }

  mdp_measure log(mdp_measure a)
  {
    return a.log();
  }

  mdp_measure pow(mdp_measure a, mdp_real b)
  {
    mdp_real mean = std::pow(a.getmean(), b);
    mdp_real error = a.getmerr() * b * std::pow(a.getmean(), b - 1);

    return mdp_measure(mean, error);
  }

  mdp_measure sin(mdp_measure a)
  {
    return a.sin();
  }

  mdp_measure cos(mdp_measure a)
  {
    return a.cos();
  }

  void print(mdp_measure a)
  {
    std::cout << a.getmean() << " (" << a.getmerr() << ")" << std::endl;
  }
} // namespace MDP

#endif /* MDP_MEASURE_ */
