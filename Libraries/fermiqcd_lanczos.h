/////////////////////////////////////////////////////////////////
/// @file fermiqcd_lanczos.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Lanczos routine
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_LANCZOS_
#define FERMIQCD_LANCZOS_

#include "mdp_complex.h"
#include "fermiqcd_gauge_field.h"
#include "fermiqcd_coefficients.h"

namespace MDP
{
  /**
   * @brief Lanczos iterative eigenvalue algorithm for lattice fermion operators.
   *
   * This class implements a single-step Lanczos iteration for computing
   * extremal eigenvalues of large sparse operators arising in lattice QCD,
   * such as Dirac-type operators.
   *
   * The algorithm builds an orthogonal Krylov basis using recursive relations:
   *
   *    A p_k = beta_{k+1} p_{k+1} + alpha_k p_k + beta_k p_{k-1}
   *
   * where A is implicitly defined through the fermion operator (mul_Q) and
   * gamma5 hermitization.
   *
   * The implementation uses a 3-vector recurrence scheme and is intended for
   * iterative use.
   *
   * @tparam fieldT Type representing lattice fermion field.
   *
   * @example
   * @verbatim
   * mdp_gauge U(lattice,nc);
   * fermi_field psi(lattice,nc);
   * coefficients coeff;
   * coeff["kappa"]=1.12;
   *
   * for(int k=0; k<100; k++)
   *   mdp << Lanczos::step(psi,U,coeff) << "\n";
   * @endverbatim
   *
   * @note This implementation uses static internal state and is not reentrant.
   */
  template <class fieldT>
  class Lanczos
  {
  private:
    Lanczos() = delete;
    ~Lanczos() = delete;
    Lanczos(const Lanczos &) = delete;
    Lanczos &operator=(const Lanczos &) = delete;

  public:
    /**
     * @brief Performs a single Lanczos iteration step.
     *
     * This function advances the Lanczos recursion by one step using a
     * 3-term recurrence relation applied to lattice fermion fields.
     *
     * The operator is implicitly defined via:
     * - mul_Q (fermion matrix application)
     * - gamma5 hermitization
     *
     * The iteration updates:
     * - p: previous Lanczos vector
     * - q: current vector
     * - r: auxiliary vector
     *
     * and computes the tridiagonal coefficients:
     * - alpha = diagonal element
     * - beta  = off-diagonal coupling
     *
     * @param psi Initial / working fermion field (modified in-place during init).
     * @param U Gauge field used in fermion operator.
     * @param coeff Coefficients defining fermion operator parameters (e.g. kappa).
     * @param force If true, reinitializes Lanczos basis.
     * @param output_check If true, prints diagnostic orthogonality checks.
     *
     * @return Complex number (alpha, beta) representing Lanczos coefficients
     *         of the tridiagonal matrix.
     *
     * @note This implementation uses static internal vectors and preserves state
     *       between calls.
     *
     * @warning Not thread-safe due to static storage of Lanczos vectors.
     *
     * @warning Numerical stability may degrade due to loss of orthogonality
     *          in long runs (no reorthogonalization is performed).
     *
     * @complexity One application of fermion operator per iteration.
     */
    static mdp_complex step(fieldT &psi,
                            gauge_field &U,
                            coefficients &coeff,
                            bool force = false,
                            bool output_check = false)
    {
      begin_function("Lanczos__step");
      static bool init = true;
      static fieldT p(psi);
      static fieldT q(psi);
      static fieldT r(psi);
      static mdp_real alpha, beta;
      static mdp_site x(psi.lattice());

      if (init || force)
      {
        mdp << "Initializing Lanczos vectors\n";

        // ///////////////////////////
        // initialize Lanczos vectors
        // ///////////////////////////

        psi.update();
        q = psi;
        mdp_real norm = std::sqrt(norm_square(q));
        q /= norm;
        p = 0.0;
        r = 0.0;
        beta = 1.0;
        init = false;
      }

      // ///////////////////////////
      // here the Lanczos algorithms
      // ///////////////////////////

      r = p;
      q /= beta;
      p = q;
      r *= -beta;
      q = r;

      p.update();
      mul_Q(r, p, U, coeff);
      multiply_by_gamma5(r, r);
      q += r;
      alpha = real_scalar_product(p, q);
      mdp_add_scaled_field(q, -alpha, p);
      beta = std::sqrt(norm_square(q));

      if (output_check)
      {
        // this prints some data for check
        static fieldT s(psi);
        mdp_real pp, qq;
        mdp_complex pq, qr, ps;
        mdp_site x(psi.lattice());

        pp = norm_square(p);
        qq = norm_square(q);
        pq = p * q;

        mul_Q(s, q, U, coeff);
        multiply_by_gamma5(s, s);
        s.update();
        qr = q * r;

        mul_Q(r, p, U, coeff);
        multiply_by_gamma5(r, r);
        r.update();
        ps = p * s;

        mdp << "|p|=" << pp << "\n";
        mdp << "|q|=" << qq << "\n";
        mdp << "&ltp|q&gt=" << pq << "\n";
        mdp << "&ltq|Q|p&gt=" << qr << "\n";
        mdp << "&ltp|Q|q&gt^*=" << ps << "\n";
        mdp << "alpha=" << alpha << "\n";
        mdp << "beta=" << beta << "\n";
      }
      end_function("Lanczos__step");
      return mdp_complex(alpha, beta);
    }
  };
} // namespace MDP

#endif /* FERMIQCD_LANCZOS_ */
