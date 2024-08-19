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
#include "mdp_prng.h"

namespace MDP
{
  /// @brief implements error propagation
  ///
  /// Example:
  /// @verbatim
  ///    mdp_measure m;
  ///    // store 10 measurements
  ///    for(int i=0; i<10; i++)
  ///       m << 3.0+mdp_random.gaussian(2.0);
  ///    m=sin(exp(m)+m);
  ///    cout << m.getmean() << "+/-" << m.getmerr() << endl;
  /// @endverbatim
  /// Assumes gaussian error propagation
  class mdp_measure
  {
  private:
    float m_mean;
    float m_error;
    int m_num;

  public:
    mdp_measure(float mean_ = 0, float error_ = 0, int num_ = 1) : m_mean(mean_), m_error(error_), m_num(num_)
    {
    }

    int getnum() const
    {
      return m_num;
    }

    float getmean() const
    {
      return m_mean;
    }

    void setmean(float mean)
    {
      m_mean = mean;
    }

    float getmerr() const
    {
      return m_error;
    }

    void setmerror(float err)
    {
      m_error = err;
    }

    void reset()
    {
      m_num = 0;
      m_mean = 0;
      m_error = 0;
    }

    void set(float x, float dx, int i = 1)
    {
      m_num = i;
      m_mean = x;
      m_error = dx;
    }

    void operator<<(float x)
    {
      float err2 = m_num * (std::pow((double)m_error, (double)2.0) + m_mean * m_mean) + std::pow((double)x, (double)2.0);
      ++m_num;
      m_mean = (m_mean * (m_num - 1) + x) / m_num;
      m_error = std::sqrt(err2 / m_num - m_mean * m_mean);
    }

    void operator>>(float &x)
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
                    std::fabs(m_mean) / std::pow((double)b.m_mean, (double)2.0) * b.m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator+(float b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean + b;
      tmp.m_error = m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator-(float b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean - b;
      tmp.m_error = m_error;
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator*(float b)
    {
      mdp_measure tmp;
      tmp.m_mean = m_mean * b;
      tmp.m_error = m_error * std::fabs(b);
      tmp.m_num = 1;
      return tmp;
    }

    mdp_measure operator/(float b)
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

  mdp_measure operator+(float a, mdp_measure b)
  {
    return b + a;
  }

  mdp_measure operator-(float a, mdp_measure b)
  {
    return -(b - a);
  }

  mdp_measure operator*(float a, mdp_measure b)
  {
    return b * a;
  }

  mdp_measure operator/(float a, mdp_measure b)
  {
    float mean = a / b.getmean();
    float error = std::fabs(a) / std::pow((double)b.getmean(), (double)2.0) * b.getmerr();
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

  mdp_measure pow(mdp_measure a, float b)
  {
    float mean = std::pow((double)a.getmean(), (double)b);
    float error = a.getmerr() * b * std::pow((double)a.getmean(), (double)b - 1);

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
