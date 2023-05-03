/*

  MDP_PVegas.h  (requires mdp.h)

  parallel implementation of the Vegas algorithm

  invented by Peter Lepage @ Cornell University

  from: vegas_int.c by Sinisa Veseli @ Fermilab, 12/1996

  rewritten in C++ by Massimo Di Pierro @ Fermilab, 4/2000

  A number of errors corrected, in particular loops now
  run from 0 to n-1, while in the original version they
  run from 1 to n, because translated badly from fortran.
  Moreover it does not need numerical recipes any more but relies
  on the Marsagla random number generator in mdp.h
  The output is much nicer now.
  Moreover all the variables and functions required are
  local or private members of VegasBase.

  Many of the original variable names have been left
  unchanged from the original version.

  Goals:
  Performs Monte Carlo integration of a user supplied
  ndim-dimensional function over a rectangular volume
  specified by { m_IntegrationLimitsMin[i], m_IntegrationLimitsMax[i] }
  The integration consists of n-iterations , each with
  approximately n-calls to the function.
  After each iteration the grid is refined; more than 5 or 10
  iterations are rarely useful. The input flag init signals
  whether this call is a new start, or a subsequent call for
  additional iterations.
  The member variable m_OutputFile (by default = stdout) is the
  (FILE*) were the basic output is directed.
  The member variable m_OutputFileAdvanced (normally = stdout) is
  the (FILE*) were the advanced output is directed.

  The main function for integration is

  double VegasClass::Integrate(int init, int ncalls, int niter)

  where function suplied in derived class is the function to integrate
        init is 0,1 or 2
        ncalls if the number of calls to function per iteration
  niter is the total number of iterations

  The output of the program is in the member variables:

  double m_Integral;
  double m_StandardDeviation;
  double m_ChiSquare;

  m_Integral is also return by Vegas.Integrate()

*/

#ifndef MDP_PVEGAS_
#define MDP_PVEGAS_

#include <algorithm>
#include <cstdio>
#include <cmath>

namespace MDP
{
  class VegasBase
  {
  protected:
    virtual double m_function(const double *) = 0;

  private:
    // constants
    static constexpr double ALPHA = 1.5;
    static constexpr int NDMX = 200;
    static constexpr int MXDIM = 10; // max grid size
    static constexpr double TINY = 1.0e-30;

    // ////////////////////////////////////////////
    // Make all variables members of the class to
    // allow restart
    // ////////////////////////////////////////////

    // variables for VegasBase:Integrate
    int m_it, m_it1, m_it2, m_mds, m_nd, m_ndo, m_ng, m_npg;
    int m_ia[MXDIM], m_kg[MXDIM];
    double m_calls, m_dv2g, m_dxg, m_ti, m_tsi, m_xjac, m_xnd;
    double m_schi, m_si, m_swgt;
    double m_d[NDMX][MXDIM], m_di[NDMX][MXDIM];
    double m_dt[MXDIM], m_dx[MXDIM];
    double m_r[NDMX], m_x[MXDIM];
    double m_xi[MXDIM][NDMX], m_xin[NDMX];
    double m_Integral;
    double m_StandardDeviation;
    enum
    {
      RelativePrecision,
      AbsolutePrecision
    } m_ConvergenceCriteria;
    double m_TargetPrecision;
    double m_ChiSquare;
    FILE *m_OutputFile;
    FILE *m_OutputFileAdvanced;
    double m_IntegrationLimitsMin[MXDIM];
    double m_IntegrationLimitsMax[MXDIM];
    int m_NumberOfDimensions;

    void PrintInputParameters(int init, int ncalls, int niterations)
    {
      if (m_OutputFile && isMainProcess())
      {
        fprintf(m_OutputFile, "======================================\n");
        fprintf(m_OutputFile, "Vegas Montecarlo Numerical Integration\n");
        fprintf(m_OutputFile, "======================================\n");
        fprintf(m_OutputFile, "Input Parameters:\n");
        fprintf(m_OutputFile, "Number of dimensions     = %i\n", m_NumberOfDimensions);
        fprintf(m_OutputFile, "Present iteration        = %i\n", m_it);
        fprintf(m_OutputFile, "Number of iterations     = %i\n", niterations);
        fprintf(m_OutputFile, "Number of function calls = %i\n", ncalls);
        fprintf(m_OutputFile, "ALPHA                    = %f\n", ALPHA);
        fprintf(m_OutputFile, "mds(?)                   = %i\n", m_mds);
        fprintf(m_OutputFile, "nd                       = %i\n\n", m_nd);
        fprintf(m_OutputFile, "Integration Limits: i\tmin(x[i])\tmax(x[i])\n");
        fprintf(m_OutputFile, "=================================================\n");
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          fprintf(m_OutputFile, "                    %i\t%f\t%f\n", j, m_IntegrationLimitsMin[j], m_IntegrationLimitsMax[j]);
        }
        fprintf(m_OutputFile, "\n");
        fflush(m_OutputFile);
      }
    }

    void PrintOutputParameters()
    {
      if (m_OutputFile && isMainProcess())
      {
        if (m_it == m_it1)
        {
          fprintf(m_OutputFile, "Iteration    Integral    StandardDeviation ChiSquare\n");
          fprintf(m_OutputFile, "====================================================\n");
        }
        fprintf(m_OutputFile, "%3d\t%14.7f\t%14.7f\t%10.3f\n", m_it, m_Integral, m_StandardDeviation, m_ChiSquare);
        if (m_it == m_it2 - 1)
        {
          fprintf(m_OutputFile, "====================================================\n");
        }
        fflush(m_OutputFile);
        PrintGrid();
      }
    }

    // /////////////////////////////////////////////////////
    // Utility routine used by integrate(), to rebin
    // a vector of densities xi into new bins defined by a
    // vector r.
    // /////////////////////////////////////////////////////
    void rebin(double rc, int nd, double r[], double xin[], double xi[])
    {
      int k = 0;
      double xo = 0.0;
      double dr = 0.0;
      double xn = 0.0;
      for (int i = 0; i < nd - 1; i++)
      {
        while (rc > dr)
        {
          dr += r[k];
          xo = xn;
          xn = xi[k];
          k++;
        }
        dr -= rc;
        xin[i] = xn - (xn - xo) * dr / r[k - 1];
      }
      for (int i = 0; i < nd - 1; i++)
      {
        xi[i] = xin[i];
      }
      xi[nd - 1] = 1.0;
    }

    int is_local(int k, int nk)
    {
      return (k % Nproc) == ME;
    }

    bool converged()
    {
      switch (m_ConvergenceCriteria)
      {
      case RelativePrecision:
        return (fabs(m_StandardDeviation / m_Integral) <= m_TargetPrecision);
      case AbsolutePrecision:
        return (fabs(m_StandardDeviation) <= m_TargetPrecision);
      default:
        return false;
      }
    }

  public:
    VegasBase()
    {
      m_OutputFile = stdout;
      m_OutputFileAdvanced = 0;
      m_ConvergenceCriteria = RelativePrecision;
      m_TargetPrecision = 1e-4;
    }

    void setDimension(int dim)
    {
      m_NumberOfDimensions = dim;
    }

    void setIntegrationLimits(int axis, int lowerBound, int upperBound)
    {
      m_IntegrationLimitsMin[axis] = lowerBound;
      m_IntegrationLimitsMax[axis] = upperBound;
    }

    void SetConvergenceCriteriaToRelativePrecision(double x = 1e-4)
    {
      m_ConvergenceCriteria = RelativePrecision;
      m_TargetPrecision = x;
    }

    void SetConvergenceCriteriaToAbsolutePrecision(double x = 1e-4)
    {
      m_ConvergenceCriteria = AbsolutePrecision;
      m_TargetPrecision = x;
    }

    void SaveGrid(const char filename[])
    {
      if (isMainProcess())
      {
        FILE *fp = fopen(filename, "w");
        fwrite(this, sizeof(VegasBase) / sizeof(char), 1, fp);
        fclose(fp);
      }
    }

    void LoadGrid(const char filename[])
    {
      FILE *fp = fopen(filename, "r");
      fread(this, sizeof(VegasBase) / sizeof(char), 1, fp);
      fclose(fp);
    }

    void PrintGrid()
    {
      if (m_OutputFileAdvanced && isMainProcess())
      {
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          fprintf(m_OutputFileAdvanced, "DATA FOR axis %2d\n", j);
          fprintf(m_OutputFileAdvanced, "Dim\ti\t x\t\tDelta_i\n");
          fprintf(m_OutputFileAdvanced, "=============================================\n");
          for (int i = 0; i < m_nd; i++)
          {
            fprintf(m_OutputFileAdvanced, "%d\t%d\t%8.5f\t%12.4f\n", j, i, m_xi[j][i], m_di[i][j]);
          }
        }
      }
      fflush(m_OutputFile);
    }

    void PrintOutput()
    {
      if (m_OutputFile && isMainProcess())
      {
        fprintf(m_OutputFile, "Integral           = %f\n", m_Integral);
        fprintf(m_OutputFile, "Standard Deviation = %f\n", m_StandardDeviation);
        fprintf(m_OutputFile, "Chi^2              = %f\n", m_ChiSquare);
        fflush(m_OutputFile);
      }
    }

    void parallel_loop(double &fb, double &f2b)
    {
      static double wgt, xn, xo, rc, f, f2;
      static double local_d[NDMX][MXDIM], local_di[NDMX][MXDIM];

      // /////////////////
      // initialization...
      // /////////////////

      fb = f2b = 0.0;
      for (int i = 0; i < m_nd; i++)
        for (int j = 0; j < m_NumberOfDimensions; j++)
          local_d[i][j] = local_di[i][j] = 0.0;

      // /////////////////
      // parallel loop
      // /////////////////
      for (int k = 0; k < m_npg; k++)
        if (is_local(k, m_npg))
        {
          wgt = m_xjac;
          for (int j = 0; j < m_NumberOfDimensions; j++)
          {
            xn = (m_kg[j] - mdp_random.plain()) * m_dxg;
            m_ia[j] = std::max(std::min((int)(xn), NDMX - 1), 0); // check NDMX
            if (m_ia[j] > 0)
            {
              xo = m_xi[j][m_ia[j]] - m_xi[j][m_ia[j] - 1];
              rc = m_xi[j][m_ia[j] - 1] + (xn - m_ia[j]) * xo;
            }
            else
            {
              xo = m_xi[j][m_ia[j]];
              rc = (xn - m_ia[j]) * xo;
            }
            m_x[j] = m_IntegrationLimitsMin[j] + rc * m_dx[j];
            wgt *= xo * m_xnd;
          }
          f = wgt * m_function(m_x);
          f2 = f * f;
          fb += f;
          f2b += f2;
          for (int j = 0; j < m_NumberOfDimensions; j++)
          {
            local_di[m_ia[j]][j] += f;
            if (m_mds >= 0)
              local_d[m_ia[j]][j] += f2;
          }
        }
      // ////////////////////////
      // Parallel sum here!
      // ////////////////////////
      mdp.add(fb);
      mdp.add(f2b);
      mdp.add((double *)local_d, m_nd * MXDIM);
      mdp.add((double *)local_di, m_nd * MXDIM);
      for (int j = 0; j < m_NumberOfDimensions; j++)
        for (int i = 0; i < m_nd; i++)
        {
          m_d[i][j] += local_d[i][j];
          m_di[i][j] += local_di[i][j];
        }
    }

    double Integrate(int init = 0,
                     int ncalls = 1000,
                     int niterations = 100)
    {
      int k;
      double wgt, f2b, fb, xo, xn, rc;

      if (init <= 0)
      {
        // ////////////////////////////////////////////
        // Normal entry. Enter here on a cold start.
        // To disable stratified sampling, i.e. to use
        // importance sampling only, change to mds = 0.
        // ////////////////////////////////////////////
        m_mds = m_ndo = 1;
        for (int j = 0; j < m_NumberOfDimensions; j++)
          m_xi[j][0] = 1.0;
      }
      if (init <= 1)
      {
        // ////////////////////////////////////////////
        // Enter here to inherit the grid from a previous
        // call, but not its answers.
        // ////////////////////////////////////////////
        m_it = 0;
        m_si = m_swgt = m_schi = 0.0;
      }
      if (init <= 2)
      {
        // ////////////////////////////////////////////
        // Enter here to inherit the previous grid
        // and its answers.
        // ////////////////////////////////////////////
        m_nd = NDMX;
        m_ng = 1;
        if (m_mds)
        { // Set up for stratification.
          m_ng = (int)std::pow(ncalls / 2.0 + 0.25, 1.0 / m_NumberOfDimensions);
          m_mds = 1;
          if ((2 * m_ng - NDMX) >= 0)
          {
            m_mds = -1;
            m_npg = m_ng / NDMX + 1;
            m_nd = m_ng / m_npg;
            m_ng = m_npg * m_nd;
          }
        }
        k = 1;
        for (int i = 0; i < m_NumberOfDimensions; i++)
        {
          k *= m_ng;
        }
        m_npg = std::max(ncalls / k, 2);
        m_calls = m_npg * k;
        m_dxg = 1.0 / m_ng;
        m_dv2g = 1;
        for (int i = 0; i < m_NumberOfDimensions; i++)
        {
          m_dv2g *= m_dxg;
        }
        m_dv2g = std::pow(m_calls * m_dv2g / m_npg, 2.0) / (m_npg - 1.0);
        m_xnd = m_nd;
        m_dxg *= m_xnd;
        m_xjac = 1.0 / m_calls;
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          m_dx[j] = m_IntegrationLimitsMax[j] - m_IntegrationLimitsMin[j];
          m_xjac *= m_dx[j];
        }
        // ////////////////////////////////////////////
        // Do binning if necessary.
        // ////////////////////////////////////////////
        if (m_nd != m_ndo)
        {
          for (int i = 0; i < m_nd; i++)
          {
            m_r[i] = 1.0;
          }
          for (int j = 0; j < m_NumberOfDimensions; j++)
          {
            rebin(m_ndo / m_xnd, m_nd, m_r, m_xin, m_xi[j]);
          }
          m_ndo = m_nd;
        }
        PrintInputParameters(init, ncalls, niterations);
      }

      // ////////////////////////////////////////////
      // Main iteration loop. Can enter here (init >= 3)
      // to do an additional niterations iterations with
      // all other parameters unchanged.
      // ////////////////////////////////////////////
      m_it1 = m_it;
      m_it2 = m_it + niterations;
      for (; m_it < m_it2; m_it++)
      {
        m_ti = m_tsi = 0.0;
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          m_kg[j] = 1;
          for (int i = 0; i < m_nd; i++)
            m_d[i][j] = m_di[i][j] = 0.0;
        }
        do
        {
          parallel_loop(fb, f2b);
          f2b = std::sqrt(f2b * m_npg);
          f2b = (f2b - fb) * (f2b + fb);
          if (f2b <= 0.0)
            f2b = TINY;
          m_ti += fb;
          m_tsi += f2b;

          if (m_mds < 0)
          { // Use stratified sampling.
            for (int j = 0; j < m_NumberOfDimensions; j++)
              m_d[m_ia[j]][j] += f2b;
          }
          for (k = m_NumberOfDimensions - 1; k >= 0; k--)
          {
            m_kg[k] %= m_ng;
            if (++m_kg[k] != 1)
              break;
          }
          if (k < 0)
            break;
        } while (true);
        // ////////////////////////////////////////////
        // Compute final results for this iteration.
        // ////////////////////////////////////////////
        m_tsi *= m_dv2g;
        wgt = 1.0 / m_tsi;

        m_si += wgt * m_ti;
        m_schi += wgt * m_ti * m_ti;
        m_swgt += wgt;
        m_Integral = m_si / m_swgt;
        if (m_it == 0)
          m_ChiSquare = 0.0;
        else
          m_ChiSquare = (m_schi - m_si * (m_Integral)) / m_it;
        if (m_ChiSquare < 0.0)
          m_ChiSquare = 0.0;
        m_StandardDeviation = std::sqrt(1.0 / m_swgt);
        m_tsi = std::sqrt(m_tsi);
        PrintOutputParameters();

        // ////////////////////////////////////////////
        // Refine the grid. The refinement is damped,
        // to avoid rapid, destabilizing changes, and also
        // compressed in range by the exponent ALPHA.
        // ////////////////////////////////////////////
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          xo = m_d[0][j];
          xn = m_d[1][j];
          m_d[0][j] = (xo + xn) / 2.0;
          m_dt[j] = m_d[0][j];
          for (int i = 1; i < m_nd - 1; i++)
          {
            rc = xo + xn;
            xo = xn;
            xn = m_d[i + 1][j];
            m_d[i][j] = (rc + xn) / 3.0;
            m_dt[j] += m_d[i][j];
          }
          m_d[m_nd - 1][j] = (xo + xn) / 2.0;
          m_dt[j] += m_d[m_nd - 1][j];
        }
        for (int j = 0; j < m_NumberOfDimensions; j++)
        {
          rc = 0.0;
          for (int i = 0; i < m_nd; i++)
          {
            if (m_d[i][j] < TINY)
              m_d[i][j] = TINY;
            m_r[i] = std::pow((1.0 - m_d[i][j] / m_dt[j]) / (std::log(m_dt[j]) - std::log(m_d[i][j])), ALPHA);
            rc += m_r[i];
          }
          rebin(rc / m_xnd, m_nd, m_r, m_xin, m_xi[j]);
        }
        if (converged())
          return m_Integral;
      }
      if (!converged() && isMainProcess())
      {
        fprintf(m_OutputFile, "Vegas failed to reach target precision.\n");
      }
      return m_Integral;
    }
  };
} // namespace MDP

#endif /* MDP_PVEGAS_ */
