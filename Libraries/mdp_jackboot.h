/////////////////////////////////////////////////////////////////
/// @file mdp_jackboot.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_jackboot
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_JACKBOOT_
#define MDP_JACKBOOT_

#include <cmath>
#include <memory>

namespace MDP
{
  /// @brief container class for jackknife and boostrap analysis
  ///
  /// Example:
  /// @verbatim
  ///    mdp_jackboot jb(10,2);
  ///    for(int k=0; k<10; k++) {
  ///       jb(k,0)=mdp_random.plain();
  ///       jb(k,1)=mdp_random.plain();
  ///    }
  ///    jb.plain(0);
  ///    cout << "mean of jb(k,0) =" << mean() << endl;
  ///    cout << "jackknife of mean jb(k,0) =" << j_err() << endl;
  ///    cout << "boostrap of mean jb(k,0) =" << b_err() << endl;
  ///    jb.plain(1);
  ///    cout << "mean of jb(k,1) =" << mean() << endl;
  ///    cout << "jackknife of mean jb(k,1) =" << j_err() << endl;
  ///    cout << "boostrap of mean jb(k,1) =" << b_err() << endl;
  /// @endverbatim
  class mdp_jackboot
  {
  private:
    int nconf;
    int narg;
    std::unique_ptr<float[]> m_data;
    int mdp_jackboot_plain_int;
    void *handle;

    void dimension(int nconf_, int narg_ = 1)
    {
      nconf = nconf_;
      narg = narg_;
      conf = 0;
      handle = nullptr;

      m_data = std::make_unique<float[]>(nconf * narg);

      for (int i = 0; i < nconf; i++)
        for (int j = 0; j < narg; j++)
          m_data[i * narg + j] = 0;
    }

    void makesample(int *p, int nboot)
    {
      for (int boot = 0; boot < nboot; boot++)
        for (int j = 0; j <= conf; j++)
          p[j + (conf + 1) * boot] = (int)((conf + 1) * mdp_random.plain());
    }

    static float mdp_jackboot_plain(float *x, void *a)
    {
      int *i_ptr = (int *)a;
      return x[*i_ptr];
    }

  public:
    float (*f)(float *, void *);
    int conf;

    /** @brief allocate container for nconf_ datasets of nargs_ numbers each
     */
    mdp_jackboot(int nconf_, int narg_ = 1)
    {
      dimension(nconf_, narg_);
    }

    virtual ~mdp_jackboot()
    {
    }

    float *address(int conf)
    {
      return m_data.get() + conf * narg;
    }

    float &operator()(int present_conf, int arg)
    {
      conf = present_conf;
      return m_data[conf * narg + arg];
    }

    float &operator()(int arg)
    {
      return m_data[conf * narg + arg];
    }

    void plain(int i)
    {
      if ((i < 0) || (i >= narg))
        error("incorrect mdp_jackboot.plain() argument");

      mdp_jackboot_plain_int = i;
      handle = (void *)&mdp_jackboot_plain_int;
      f = mdp_jackboot_plain;
    }

    float mean()
    {
      std::unique_ptr<float[]> average = std::make_unique<float[]>(narg);

      if (conf == 0)
        return (*f)(&m_data[0], handle);

      for (int i = 0; i < narg; i++)
      {
        average[i] = 0;
        for (int j = 0; j <= conf; j++)
          average[i] += m_data[j * narg + i] / (conf + 1);
      }

      return (*f)(average.get(), handle);
    }

    float j_err()
    {
      std::unique_ptr<float[]> average = std::make_unique<float[]>((conf + 1) * narg);

      if (conf < 2)
        return 0;

      for (int i = 0; i < narg; i++)
        for (int j = 0; j <= conf; j++)
        {
          average[j * narg + i] = 0;
          for (int k = 0; k <= conf; k++)
            // typo noticed by Chris Schroeder was fixed!
            if (j != k)
              average[j * narg + i] +=
                  (m_data[k * narg + i] / (nconf - 1));
        }

      float mean0 = mean();
      float stddev = 0;

      for (int j = 0; j <= conf; j++)
        stddev += (float)std::pow(
            (double)(*f)(&(average[j * narg]), handle) - mean0, 2.0);

      return (float)std::sqrt(stddev * conf / (conf + 1));
    }

    float b_err(int nboot = 100)
    {
      float vx, vy, tmp;
      std::unique_ptr<float[]> average = std::make_unique<float[]>(nboot * narg);
      std::unique_ptr<int[]> p = std::make_unique<int[]>((conf + 1) * nboot);

      makesample(p.get(), nboot);

      if (conf == 0)
        return 0;
      for (int i = 0; i < narg; i++)
        for (int boot = 0; boot < nboot; boot++)
        {
          average[boot * narg + i] = 0;
          for (int j = 0; j <= conf; j++)
            average[boot * narg + i] += m_data[p[j + (conf + 1) * boot] * narg + i];
          average[boot * narg + i] /= conf + 1;
        }

      for (int x = 1; x < nboot; x++)
      {
        vx = (*f)(&(average[x * narg]), handle);
        for (int y = x - 1; y >= 0; y--)
        {
          vy = (*f)(&(average[y * narg]), handle);
          if (vy > vx)
            for (int i = 0; i < narg; i++)
            {
              tmp = average[(y + 1) * narg + i];
              average[(y + 1) * narg + i] = average[y * narg + i];
              average[y * narg + i] = tmp;
            }
          else
            y = -1;
        }
      }

      vx = (*f)(&(average[((int)((float)nboot * 0.16)) * narg]), handle);
      vy = (*f)(&(average[((int)((float)nboot * 0.84)) * narg]), handle);

      return (float)(vy - vx) / 2.0;
    }
  };
} // namespace MDP

#endif /* MDP_JACKBOOT_ */
