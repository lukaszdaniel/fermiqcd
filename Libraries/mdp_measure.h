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
  ///    cout << m.getmean() << "+/-" << m.geterr() << endl;
  /// @endverbatim
  /// Assumes gaussian error propagation
  class mdp_measure
  {
  public:
    int num;
    float mean;
    float error;

    int getnum() const
    {
      return num;
    }

    float getmean() const
    {
      return mean;
    }

    float getmerr() const
    {
      return error;
    }

    mdp_measure()
    {
      num = 0;
      mean = 0;
      error = 0;
    }

    mdp_measure(float mean_, float error_, int num_ = 1)
    {
      num = num_;
      mean = mean_;
      error = error_;
    }

    void reset()
    {
      num = 0;
      mean = 0;
      error = 0;
    }

    void set(float x, float dx, int i = 1)
    {
      num = i;
      mean = x;
      error = dx;
    }

    void operator<<(float x)
    {
      float err2 = num * (std::pow((double)error, (double)2.0) + mean * mean) + std::pow((double)x, (double)2.0);
      ++num;
      mean = (mean * (num - 1) + x) / num;
      error = std::sqrt(err2 / num - mean * mean);
    }

    void operator>>(float &x)
    {
      x = mean + error * Random.gaussian();
    }

    friend mdp_measure operator+(mdp_measure a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean + b.mean;
      tmp.error = a.error + b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator-(mdp_measure a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean - b.mean;
      tmp.error = a.error + b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator*(mdp_measure a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean * b.mean;
      tmp.error = std::fabs(a.mean) * b.error + a.error * std::fabs(b.mean);
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator/(mdp_measure a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean / b.mean;
      tmp.error = a.error / std::fabs(b.mean) +
                  std::fabs(a.mean) / std::pow((double)b.mean, (double)2.0) * b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator+(float a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a + b.mean;
      tmp.error = b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator-(float a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a - b.mean;
      tmp.error = b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator*(float a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a * b.mean;
      tmp.error = std::fabs(a) * b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator/(float a, mdp_measure b)
    {
      mdp_measure tmp;
      tmp.mean = a / b.mean;
      tmp.error = std::fabs(a) / std::pow((double)b.mean, (double)2.0) * b.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator+(mdp_measure a, float b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean + b;
      tmp.error = a.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator-(mdp_measure a, float b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean - b;
      tmp.error = a.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator*(mdp_measure a, float b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean * b;
      tmp.error = a.error * std::fabs(b);
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure operator/(mdp_measure a, float b)
    {
      mdp_measure tmp;
      tmp.mean = a.mean / b;
      tmp.error = a.error / std::fabs(b);
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure exp(mdp_measure a)
    {
      mdp_measure tmp;
      tmp.mean = std::exp(a.mean);
      tmp.error = std::exp(a.mean) * a.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure log(mdp_measure a)
    {
      mdp_measure tmp;
      tmp.mean = std::log(a.mean);
      tmp.error = a.error / fabs(a.mean);
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure pow(mdp_measure a, float b)
    {
      mdp_measure tmp;
      tmp.mean = std::pow((double)a.mean, (double)b);
      tmp.error = a.error * b * std::pow((double)a.mean, (double)b - 1);
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure sin(mdp_measure a)
    {
      mdp_measure tmp;
      tmp.mean = std::sin(a.mean);
      tmp.error = std::fabs(std::cos(a.mean)) * a.error;
      tmp.num = 1;
      return tmp;
    }

    friend mdp_measure cos(mdp_measure a)
    {
      mdp_measure tmp;
      tmp.mean = std::cos(a.mean);
      tmp.error = std::fabs(std::sin(a.mean)) * a.error;
      tmp.num = 1;
      return tmp;
    }

    friend void print(mdp_measure a)
    {
      printf("%f (%f)\n", a.mean, a.error);
    }
  };
} // namespace MDP

#endif /* MDP_MEASURE_ */
