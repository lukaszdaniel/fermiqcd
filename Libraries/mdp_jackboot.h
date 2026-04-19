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
#include "mdp_global_vars.h"
#include "mdp_macros.h"
#include "mdp_random.h"

namespace MDP
{
  /// @brief container class for jackknife and boostrap analysis
  ///
  /// Example:
  /// @verbatim
  ///    mdp_jackboot jb(10,2);
  ///    for(mdp_uint k=0; k<10; k++) {
  ///       jb(k,0)=mdp_global_random.plain();
  ///       jb(k,1)=mdp_global_random.plain();
  ///    }
  ///    jb.plain(0);
  ///    std::cout << "mean of jb(k,0) =" << mean() << std::endl;
  ///    std::cout << "jackknife of mean jb(k,0) =" << j_err() << std::endl;
  ///    std::cout << "boostrap of mean jb(k,0) =" << b_err() << std::endl;
  ///    jb.plain(1);
  ///    std::cout << "mean of jb(k,1) =" << mean() << std::endl;
  ///    std::cout << "jackknife of mean jb(k,1) =" << j_err() << std::endl;
  ///    std::cout << "boostrap of mean jb(k,1) =" << b_err() << std::endl;
  /// @endverbatim
  class mdp_jackboot
  {
  private:
    mdp_uint m_nconf;
    mdp_uint m_narg;
    mdp_uint m_conf;
    std::unique_ptr<mdp_real[]> m_data;
    mdp_uint m_mdp_jackboot_plain_int;
    const void *m_handle;

    /** @brief Initialize jackboot dataset
     *
     * @param nconf_ Number of datasets
     * @param narg_ Number of variables
     */
    void dimension(mdp_uint nconf_, mdp_uint narg_ = 1)
    {
      m_nconf = nconf_;
      m_narg = narg_;
      m_conf = 0;
      m_handle = nullptr;

      m_data = std::make_unique<mdp_real[]>(m_nconf * m_narg);

      for (mdp_uint i = 0; i < m_nconf; i++)
        for (mdp_uint j = 0; j < m_narg; j++)
          m_data[i * m_narg + j] = 0;
    }

    void makesample(mdp_uint p[], mdp_uint nboot)
    {
      for (mdp_uint boot = 0; boot < nboot; boot++)
        for (mdp_uint j = 0; j <= m_conf; j++)
          p[j + (m_conf + 1) * boot] = (mdp_uint)((m_conf + 1) * mdp_global_random.plain());
    }

    static mdp_real mdp_jackboot_plain(const mdp_real x[], const void *a)
    {
      const mdp_uint *i_ptr = static_cast<const mdp_uint *>(a);
      return x[*i_ptr];
    }

  public:
    mdp_real (*f)(const mdp_real *, const void *);

    /** @brief allocate container for nconf_ datasets of nargs_ numbers each
     *
     * @param nconf_ Number of datasets
     * @param narg_ Number of variables
     */
    mdp_jackboot(mdp_uint nconf_, mdp_uint narg_ = 1)
    {
      dimension(nconf_, narg_);
    }

    ~mdp_jackboot()
    {
    }

    mdp_real *address(mdp_uint conf)
    {
      return m_data.get() + conf * m_narg;
    }

    mdp_real &operator()(mdp_uint present_conf, mdp_uint arg)
    {
      m_conf = present_conf;
      return m_data[m_conf * m_narg + arg];
    }

    mdp_real &operator()(mdp_uint arg)
    {
      return m_data[m_conf * m_narg + arg];
    }

    void set_conf(mdp_uint conf)
    {
      m_conf = conf;
    }

    void plain(mdp_uint i)
    {
      m_mdp_jackboot_plain_int = i;
      m_handle = (const void *)&m_mdp_jackboot_plain_int;
      f = mdp_jackboot_plain;
    }

    /** @brief Calculate sample mean
     */
    mdp_real mean()
    {
      std::unique_ptr<mdp_real[]> average = std::make_unique<mdp_real[]>(m_narg);

      if (m_conf == 0)
        return f(&m_data[0], m_handle);

      for (mdp_uint i = 0; i < m_narg; i++)
      {
        average[i] = 0;
        for (mdp_uint j = 0; j <= m_conf; j++)
          average[i] += m_data[j * m_narg + i] / (m_conf + 1);
      }

      return f(average.get(), m_handle);
    }

    /** @brief Calculate Jacknife error
     */
    mdp_real j_err()
    {
      std::unique_ptr<mdp_real[]> average = std::make_unique<mdp_real[]>((m_conf + 1) * m_narg);

      if (m_conf < 2)
        return 0;

      for (mdp_uint i = 0; i < m_narg; i++)
        for (mdp_uint j = 0; j <= m_conf; j++)
        {
          average[j * m_narg + i] = 0;
          for (mdp_uint k = 0; k <= m_conf; k++)
            // typo noticed by Chris Schroeder was fixed!
            if (j != k)
              average[j * m_narg + i] +=
                  (m_data[k * m_narg + i] / (m_nconf - 1));
        }

      mdp_real mean0 = mean();
      mdp_real stddev = 0;

      for (mdp_uint j = 0; j <= m_conf; j++)
        stddev += std::pow(f(&(average[j * m_narg]), m_handle) - mean0, 2.0);

      return std::sqrt(stddev * m_conf / (m_conf + 1));
    }

    /** @brief Calculate Bootstrap error
     *
     * @param nboot Number of boostrap iterations
     */
    mdp_real b_err(mdp_uint nboot = 100)
    {
      mdp_real vx, vy;
      std::unique_ptr<mdp_real[]> average = std::make_unique<mdp_real[]>(nboot * m_narg);
      std::unique_ptr<mdp_uint[]> p = std::make_unique<mdp_uint[]>((m_conf + 1) * nboot);

      makesample(p.get(), nboot);

      if (m_conf == 0)
        return 0;

      for (mdp_uint i = 0; i < m_narg; i++)
        for (mdp_uint boot = 0; boot < nboot; boot++)
        {
          average[boot * m_narg + i] = 0;
          for (mdp_uint j = 0; j <= m_conf; j++)
            average[boot * m_narg + i] += m_data[p[j + (m_conf + 1) * boot] * m_narg + i];
          average[boot * m_narg + i] /= m_conf + 1;
        }

      for (mdp_uint x = 1; x < nboot; x++)
      {
        vx = f(&(average[x * m_narg]), m_handle);
        for (mdp_uint y = x; y > 0; y--)
        {
          vy = f(&(average[(y - 1) * m_narg]), m_handle);
          if (vy > vx)
            for (mdp_uint i = 0; i < m_narg; i++)
            {
              std::swap(average[y * m_narg + i], average[(y - 1) * m_narg + i]);
              // mdp_real tmp = average[y * m_narg + i];
              // average[y * m_narg + i] = average[(y - 1) * m_narg + i];
              // average[(y - 1) * m_narg + i] = tmp;
            }
          else
            y = -1;
        }
      }

      vx = f(&(average[(mdp_uint(nboot * 0.16)) * m_narg]), m_handle);
      vy = f(&(average[(mdp_uint(nboot * 0.84)) * m_narg]), m_handle);

      return (vy - vx) / 2.0;
    }
  };
} // namespace MDP

#endif /* MDP_JACKBOOT_ */
