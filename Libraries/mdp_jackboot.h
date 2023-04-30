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
    int m_nconf;
    int m_narg;
    int m_conf;
    std::unique_ptr<float[]> m_data;
    int m_mdp_jackboot_plain_int;
    const void *m_handle;

    /** @brief Initialize jackboot dataset
     *
     * @param nconf_ Number of datasets
     * @param narg_ Number of variables
     */
    void dimension(int nconf_, int narg_ = 1)
    {
      m_nconf = nconf_;
      m_narg = narg_;
      m_conf = 0;
      m_handle = nullptr;

      m_data = std::make_unique<float[]>(m_nconf * m_narg);

      for (int i = 0; i < m_nconf; i++)
        for (int j = 0; j < m_narg; j++)
          m_data[i * m_narg + j] = 0;
    }

    void makesample(int *p, int nboot)
    {
      for (int boot = 0; boot < nboot; boot++)
        for (int j = 0; j <= m_conf; j++)
          p[j + (m_conf + 1) * boot] = (int)((m_conf + 1) * mdp_random.plain());
    }

    static float mdp_jackboot_plain(const float *x, const void *a)
    {
      const int *i_ptr = static_cast<const int *>(a);
      return x[*i_ptr];
    }

  public:
    float (*f)(const float *, const void *);

    /** @brief allocate container for nconf_ datasets of nargs_ numbers each
     *
     * @param nconf_ Number of datasets
     * @param narg_ Number of variables
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
      return m_data.get() + conf * m_narg;
    }

    float &operator()(int present_conf, int arg)
    {
      m_conf = present_conf;
      return m_data[m_conf * m_narg + arg];
    }

    float &operator()(int arg)
    {
      return m_data[m_conf * m_narg + arg];
    }

    void set_conf(int conf)
    {
      m_conf = conf;
    }

    void plain(int i)
    {
      if ((i < 0) || (i >= m_narg))
        error("incorrect mdp_jackboot.plain() argument");

      m_mdp_jackboot_plain_int = i;
      m_handle = (const void *)&m_mdp_jackboot_plain_int;
      f = mdp_jackboot_plain;
    }

    /** @brief Calculate sample mean
     */
    float mean()
    {
      std::unique_ptr<float[]> average = std::make_unique<float[]>(m_narg);

      if (m_conf == 0)
        return (*f)(&m_data[0], m_handle);

      for (int i = 0; i < m_narg; i++)
      {
        average[i] = 0;
        for (int j = 0; j <= m_conf; j++)
          average[i] += m_data[j * m_narg + i] / (m_conf + 1);
      }

      return (*f)(average.get(), m_handle);
    }

    /** @brief Calculate Jacknife error
     */
    float j_err()
    {
      std::unique_ptr<float[]> average = std::make_unique<float[]>((m_conf + 1) * m_narg);

      if (m_conf < 2)
        return 0;

      for (int i = 0; i < m_narg; i++)
        for (int j = 0; j <= m_conf; j++)
        {
          average[j * m_narg + i] = 0;
          for (int k = 0; k <= m_conf; k++)
            // typo noticed by Chris Schroeder was fixed!
            if (j != k)
              average[j * m_narg + i] +=
                  (m_data[k * m_narg + i] / (m_nconf - 1));
        }

      float mean0 = mean();
      float stddev = 0;

      for (int j = 0; j <= m_conf; j++)
        stddev += (float)std::pow(
            (double)(*f)(&(average[j * m_narg]), m_handle) - mean0, 2.0);

      return (float)std::sqrt(stddev * m_conf / (m_conf + 1));
    }

    /** @brief Calculate Bootstrap error
     *
     * @param nboot Number of boostrap iterations
     */
    float b_err(int nboot = 100)
    {
      float vx, vy, tmp;
      std::unique_ptr<float[]> average = std::make_unique<float[]>(nboot * m_narg);
      std::unique_ptr<int[]> p = std::make_unique<int[]>((m_conf + 1) * nboot);

      makesample(p.get(), nboot);

      if (m_conf == 0)
        return 0;

      for (int i = 0; i < m_narg; i++)
        for (int boot = 0; boot < nboot; boot++)
        {
          average[boot * m_narg + i] = 0;
          for (int j = 0; j <= m_conf; j++)
            average[boot * m_narg + i] += m_data[p[j + (m_conf + 1) * boot] * m_narg + i];
          average[boot * m_narg + i] /= m_conf + 1;
        }

      for (int x = 1; x < nboot; x++)
      {
        vx = (*f)(&(average[x * m_narg]), m_handle);
        for (int y = x - 1; y >= 0; y--)
        {
          vy = (*f)(&(average[y * m_narg]), m_handle);
          if (vy > vx)
            for (int i = 0; i < m_narg; i++)
            {
              tmp = average[(y + 1) * m_narg + i];
              average[(y + 1) * m_narg + i] = average[y * m_narg + i];
              average[y * m_narg + i] = tmp;
            }
          else
            y = -1;
        }
      }

      vx = (*f)(&(average[((int)((float)nboot * 0.16)) * m_narg]), m_handle);
      vy = (*f)(&(average[((int)((float)nboot * 0.84)) * m_narg]), m_handle);

      return (float)(vy - vx) / 2.0;
    }
  };
} // namespace MDP

#endif /* MDP_JACKBOOT_ */
