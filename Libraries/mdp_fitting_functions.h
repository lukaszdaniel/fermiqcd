/////////////////////////////////////////////////////////////////
/// @file mdp_fitting_functions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains fitting functions
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_FITTING_FUNCTIONS_
#define MDP_FITTING_FUNCTIONS_

#include <memory>
#include "mdp_measure.h"

namespace MDP
{
  /** @brief Fits y[i], x[i] for i0<=i<in with y=a[0]*x+a[1]
   *  Weighted least squares linear regression
   */
  void linear_fit(mdp_real *x, mdp_measure *y, mdp_int i0, mdp_int in, mdp_measure *a)
  {
    mdp_int i;
    mdp_real S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (i = i0; i < in; i++)
    {
      mdp_real dy2 = std::pow(y[i].getmerr(), 2.0);
      if (dy2 <= 0.0)
        dy2 = std::pow(y[i].getmean() * 1e-5, 2.0);

      S = S + 1.0 / dy2;
      Sx = Sx + 1.0 / dy2 * x[i];
      Sy = Sy + y[i].getmean() / dy2;
      Sxx = Sxx + x[i] * x[i] / dy2;
      Sxy = Sxy + y[i].getmean() / dy2 * x[i];
    }

    const mdp_real det = S * Sxx - Sx * Sx;

    // this is m
    a[0].setmean((S * Sxy - Sx * Sy) / det);
    a[0].setmerror(std::sqrt(S / det));
    // this is q
    a[1].setmean((Sxx * Sy - Sx * Sxy) / det);
    a[1].setmerror(std::sqrt(Sxx / det));
  }

  /** @brief Fitting function
   *  finds x=xmin that minimizes
   *  (*fp)(&x,1,dummy)
   *   requires that:
   *  (*fp)(&ax) > (*fp)(&bx) && (*fp)(&cx) > (*fp)(&bx)
   */
  mdp_real golden_rule(mdp_real (*fp)(mdp_real *, mdp_int, void *), mdp_real &xmin,
                    mdp_real ax, mdp_real bx, mdp_real cx,
                    mdp_real tol = 0.001, mdp_int niter = 100, void *dummy = nullptr)
  {
    static constexpr mdp_real CGOLD = 0.3819660;

    mdp_real a = (ax < cx ? ax : cx);
    mdp_real b = (ax > cx ? ax : cx);
    mdp_real x = bx, w = bx, v = bx;

    mdp_real fx = fp(&x, 1, dummy);
    mdp_real fw = fx, fv = fx;
    mdp_real e = 0.0;

    for (mdp_int iter = 0; iter < niter; iter++)
    {
      const mdp_real xm = 0.5 * (a + b);
      const mdp_real tol1 = tol * fabs(x) + mdp_precision;
      const mdp_real tol2 = 2.0 * tol1;

      if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
      {
        xmin = x;
        return fx;
      }

      mdp_real d = 0.0;

      if (fabs(e) > tol1)
      {
        // Parabolic interpolation
        const mdp_real r = (x - w) * (fx - fv);
        const mdp_real q = (x - v) * (fx - fw);
        mdp_real p = (x - v) * q - (x - w) * r;
        mdp_real q_adj = 2.0 * (q - r);

        if (q_adj > 0.0)
          p = -p;
        q_adj = fabs(q_adj);

        const mdp_real etemp = e;
        e = d;

        if (fabs(p) < fabs(0.5 * q_adj * etemp) && p > q_adj * (a - x) && p < q_adj * (b - x))
        {
          d = p / q_adj;
          const mdp_real u = x + d;
          if ((u - a) < tol2 || (b - u) < tol2)
          {
            d = fabs(tol1);
            if (xm - x < 0.0)
              d = -d;
          }
        }
        else
        {
          d = CGOLD * (e = ((x >= xm) ? (a - x) : (b - x)));
        }
      }
      else
      {
        d = CGOLD * (e = (x >= xm ? (a - x) : (b - x)));
      }

      mdp_real u = 0.0;
      if (fabs(d) >= tol1)
        u = x + d;
      else
        u = x + (d > 0 ? fabs(tol1) : -fabs(tol1));

      mdp_real fu = fp(&u, 1, dummy);

      if (fu <= fx)
      {
        if (u > x)
          a = x;
        else
          b = x;
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      }
      else
      {
        if (u < x)
          a = u;
        else
          b = u;
        if (fu < fw || w == x)
        {
          v = w;
          w = u;
          fv = fw;
          fw = fu;
        }
        else if (fu <= fv || v == x || v == w)
        {
          v = u;
          fv = fu;
        }
      }
    }

    xmin = x;
    return fx;
  }

  using BLM_function = mdp_real (*)(mdp_real, mdp_real *, mdp_int, void *);

  /** @brief This function is used by the BayesianLevenbergMarquardt
   * It computes the chi_square (including the Bayesian term)
   * and fills alpha and beta
   *
   * \f$ \alpha(j,k)=\sum_i (Dy(x[i],a)/Da[j])*(Dy(x[i],a)/Da[k])/dy[i]^2\f$
   *
   * \f$ \beta(j)=sum_i     (y[i]-y(x[i],a))*(dy(x[i],a)/da[j])/dy[i]^2\f$
   *
   * \f$ \chi_{square}=\sum_i (y[i]-y(x[i],a))*(y[i]-y(x[i],a))/dy[i]^2
   *           +\sum_{j,k} (a[j]-a0[j])*(a[k]-a0[k])*\sigma(j,k)\f$
   *
   * This function take into account multipliticty factors
   * y[i].num, i.e. the numbers of measures used to determine y[i].mean
   * This is used as a weight factor!
   */
  mdp_real BLMaux(mdp_real *x, mdp_measure *y,
               mdp_int i_min, mdp_int i_max,
               mdp_real *a, mdp_real *a0, mdp_matrix &sigma, int ma,
               mdp_matrix &alpha,
               mdp_matrix &beta,
               BLM_function func,
               mdp_real h,
               void *junk)
  {
    std::unique_ptr<mdp_real[]> dyda = std::make_unique<mdp_real[]>(ma);

    for (mdp_int i = 0; i < ma; i++)
    {
      for (mdp_int j = 0; j < ma; j++)
      {
        alpha(i, j) = 0.0;
      }
      beta(i, 0) = 0.0;
    }

    mdp_real chi_square = 0.0;

    for (mdp_int i = i_min; i < i_max; i++)
    {
      const mdp_real ymod = func(x[i], a, ma, junk);

      // Numeric derivatives
      for (mdp_int j = 0; j < ma; j++)
      {
        a[j] += h;
        dyda[j] = func(x[i], a, ma, junk);
        a[j] -= 2 * h;
        dyda[j] -= func(x[i], a, ma, junk);
        a[j] += h;
        dyda[j] /= (2.0 * h);
      }

      const mdp_real sig2i = 1.0 / std::pow(y[i].getmerr(), 2.0);
      const mdp_real dy = y[i].getmean() - ymod;
      const mdp_real weight = y[i].getnum();

      for (mdp_int j = 0; j < ma; j++)
      {
        const mdp_real wt = dyda[j] * sig2i;
        for (mdp_int k = 0; k < ma; k++)
        {
          alpha(j, k) += dyda[k] * wt * weight;
        }
        beta(j, 0) += dy * wt * weight;
      }

      chi_square += dy * dy * sig2i;
    }

    for (mdp_int i = 0; i < ma; i++)
    {
      for (mdp_int j = 0; j < ma; j++)
      {
        chi_square += (a0[i] - a[i]) * sigma(i, j).real() * (a0[j] - a[j]);
      }
    }

    return chi_square;
  }

  /** @brief This implements the BayesianLevenbergMarquardt
   *
   * It uses mdp_matrix.
   * Arguments are:
   *
   * x[i]          : an array of mdp_real
   * y[i]          : an array of Measures
   * i_min, i_max  : range to be used in the fit
   *                 points within the range that have y[i].num=0
   *                 are ignored
   * a[i], ma      : vector of paramters for the fit and number of parameters
   *                 they are all used in the fit
   *                 the initial values are used as preons
   * covar(i,j)    : covariance matrix for the preons
   * func(x,a,ma,junk) : the function to be used in the fit
   * h             : a mdp_real used to evaluate derivatives
   * nmax          : max number of iterations
   * junk          : junk to be passed to func
   *
   * Return the Bayesian ChiSquare. To obtain the correct chi_square
   *  rerun it with same fitting values and nmax=1;
   */
  mdp_real BayesianLevenbergMarquardt(mdp_real *x, mdp_measure *y,
                                   mdp_int i_min, mdp_int i_max,
                                   mdp_real *a, int ma,
                                   mdp_matrix &covar,
                                   BLM_function func,
                                   mdp_real h = 0.001,
                                   mdp_int nmax = 1000,
                                   void *junk = nullptr)
  {

    mdp_real lambda = 0.1;
    constexpr mdp_real scale = 2.0;

    mdp_matrix alpha(ma, ma);
    std::unique_ptr<mdp_real[]> atry = std::make_unique<mdp_real[]>(ma);
    std::unique_ptr<mdp_real[]> a0 = std::make_unique<mdp_real[]>(ma);
    mdp_matrix sigma(ma, ma);
    mdp_matrix beta(ma, 1);
    mdp_matrix oneda(ma, 1);

    if (max(covar) <= 1e-6)
    {
      for (mdp_int i = 0; i < ma; i++)
      {
        sigma(i, i) = 0;
      }
    }
    else
    {
      sigma = inv(covar);
    }

    for (mdp_int i = 0; i < ma; i++)
    {
      atry[i] = a[i];
      a0[i] = a[i];
    }

    mdp_real chi_square = BLMaux(x, y, i_min, i_max, a, a0.get(), sigma,
                              ma, alpha, beta, func, h, junk);
    mdp_real old_chi_square = chi_square;

    for (mdp_int n = 0; n < nmax; n++)
    {
      for (mdp_int i = 0; i < ma; i++)
      {
        for (mdp_int j = 0; j < ma; j++)
        {
          covar(i, j) = alpha(i, j);
        }
        covar(i, i) *= (1.0 + lambda);
        oneda(i, 0) = beta(i, 0);
      }

      // this is gaussj
      covar = inv(covar);
      oneda = covar * oneda;

      for (mdp_int i = 0; i < ma; i++)
      {
        atry[i] = a[i] + real(oneda(i, 0));
      }

      chi_square = BLMaux(x, y, i_min, i_max, atry.get(), a0.get(), sigma,
                          ma, covar, oneda, func, h, junk);

      if (fabs(chi_square - old_chi_square) < 1e-4 || lambda == 0)
      {
        covar = inv(covar);
        return chi_square;
      }

      if (chi_square <= old_chi_square)
      {
        lambda /= scale;
        old_chi_square = chi_square;

        for (mdp_int i = 0; i < ma; i++)
        {
          for (mdp_int j = 0; j < ma; j++)
          {
            alpha(i, j) = covar(i, j);
          }

          beta(i, 0) = oneda(i, 0);
        }

        for (mdp_int j = 0; j < ma; j++)
        {
          a[j] = atry[j];
        }
      }
      else
      {
        lambda *= scale;
        chi_square = old_chi_square;
      }
    }

    return chi_square;
  }
} // namespace MDP

#endif /* MDP_FITTING_FUNCTIONS_ */
