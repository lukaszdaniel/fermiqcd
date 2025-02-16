/////////////////////////////////////////////////////////////////
/// @file fermiqcd_ffts.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Discrete Fourier stransform (not FFT quote yet)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_FFTS_
#define FERMIQCD_FFTS_

#include "mdp_global_vars.h"
#include "mdp_field.h"
#include "fermiqcd_fermi_field.h"

namespace MDP
{
  mdp_int i2pow(mdp_int n)
  {
    return 0x0001 << n;
  }

  /** @brief one dimensional Fourier transform
   *
   * here n is not the size of f but size of f is i2pow(n)
   */
  void dft(mdp_complex *fft_f, mdp_complex *f, mdp_int n, double sign,
           mdp_int offset = 0, mdp_int coeff = 1)
  {
    mdp_complex phase = exp(2.0 * Pi * I * (1.0 * sign / n));
    for (int i = 0; i < n; i++)
    {
      fft_f[offset + coeff * i] = 0;
      for (int j = 0; j < n; j++)
        fft_f[offset + coeff * i] += f[offset + coeff * j] * pow(phase, i * j);
      fft_f[offset + coeff * i] /= std::sqrt(n);
    }
  }

  void fft(mdp_complex *fft_f, mdp_complex *f, mdp_int n, double sign,
           mdp_int offset = 0, mdp_int coeff = 1)
  {
    mdp_int N = i2pow(n); // Size of the transform (2^n)
    if (sign != 0.0)
    {

      // Copy data from f to fft_f
      for (mdp_int i = 0; i < N; i++)
        fft_f[offset + coeff * i] = f[offset + coeff * i];

      // Bit reversal step (important for FFT)
      mdp_int j = 0;
      for (mdp_int i = 1; i < N; i++)
      {
        mdp_int bit = N >> 1;
        for (; j >= bit; bit >>= 1)
          j -= bit;
        j += bit;

        if (i < j)
          std::swap(fft_f[offset + coeff * i], fft_f[offset + coeff * j]);
      }

      // Cooley-Tukey FFT
      for (mdp_int s = 1; s <= n; s++)
      {
        mdp_int m = std::pow(2, s); // Size of the current stage
        mdp_int half_m = m / 2;
        mdp_complex omega_m = exp(mdp_complex(0, sign * 2.0 * Pi / m)); // Twiddle factor

        for (mdp_int k = 0; k < N; k += m)
        {
          mdp_complex omega = 1.0;
          for (mdp_int j = 0; j < half_m; j++)
          {
            mdp_complex t = omega * fft_f[offset + coeff * (k + j + half_m)];
            mdp_complex u = fft_f[offset + coeff * (k + j)];
            fft_f[offset + coeff * (k + j)] = u + t;
            fft_f[offset + coeff * (k + j + half_m)] = u - t;

            omega *= omega_m; // Update the omega value for the next iteration
          }
        }
      }

      // Normalize if needed, typically for FFT we divide by sqrt(N) if not an inverse FFT
      for (mdp_int i = 0; i < N; i++)
      {
        fft_f[offset + coeff * i] /= std::sqrt(N);
      }
    }
    else
    {
      for (mdp_int a = 0; a < N; a++)
        fft_f[offset + coeff * a] = f[offset + coeff * a];
    }
  }

  void fermi_field_fft(int t,
                       fermi_field &psi_out,
                       fermi_field &psi_in,
                       int sign)
  {

    if (psi_in.lattice().n_dimensions() != 4)
      error("fft3D requires TxXxXxX");

    mdp_int size = psi_in.lattice().size(1);
    if (psi_in.lattice().size(2) > size)
      size = psi_in.lattice().size(2);
    if (psi_in.lattice().size(3) > size)
      size = psi_in.lattice().size(3);

    std::unique_ptr<mdp_complex[]> v = std::make_unique<mdp_complex[]>(size);
    std::unique_ptr<mdp_complex[]> u = std::make_unique<mdp_complex[]>(size);

    mdp_site x(psi_in.lattice());

    forallsites(x)
    {
      if (x(0) == t)
        psi_out(x) = psi_in(x);
    }

    for (mdp_int spin = 0; spin < psi_out.nspin(); spin++)
      for (mdp_int color = 0; color < psi_out.nc(); color++)
      {
        for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (mdp_int i = 0; i < psi_out.lattice().size(1); i++)
            {
              x.set(t, i, x2, x3);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(1), sign);
            for (mdp_int i = 0; i < psi_out.lattice().size(1); i++)
            {
              x.set(t, i, x2, x3);
              psi_out(x, spin, color) = u[i];
            }
          }
        for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (mdp_int i = 0; i < psi_out.lattice().size(2); i++)
            {
              x.set(t, x1, i, x3);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(2), sign);
            for (mdp_int i = 0; i < psi_out.lattice().size(2); i++)
            {
              x.set(t, x1, i, x3);
              psi_out(x, spin, color) = u[i];
            }
          }

        for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          {
            for (mdp_int i = 0; i < psi_out.lattice().size(3); i++)
            {
              x.set(t, x1, x2, i);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(3), sign);
            for (mdp_int i = 0; i < psi_out.lattice().size(3); i++)
            {
              x.set(t, x1, x2, i);
              psi_out(x, spin, color) = u[i];
            }
          }
      }
  }

  void fermi_field_fft_t(fermi_field &psi_out,
                         fermi_field &psi_in,
                         int sign)
  {

    if (psi_in.lattice().n_dimensions() != 4)
      error("fft3D requires TxXxXxX");

    mdp_int size = psi_in.lattice().size(0);
    if (psi_in.lattice().size(2) > size)
      size = psi_in.lattice().size(2);
    if (psi_in.lattice().size(3) > size)
      size = psi_in.lattice().size(3);

    std::unique_ptr<mdp_complex[]> v = std::make_unique<mdp_complex[]>(size);
    std::unique_ptr<mdp_complex[]> u = std::make_unique<mdp_complex[]>(size);

    mdp_site x(psi_in.lattice());

    forallsites(x)
    {
      psi_out(x) = psi_in(x);
    }

    for (mdp_int spin = 0; spin < psi_out.nspin(); spin++)
      for (mdp_int color = 0; color < psi_out.nc(); color++)
      {
        for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
            for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
            {
              for (mdp_int i = 0; i < psi_out.lattice().size(0); i++)
              {
                x.set(i, x1, x2, x3);
                v[i] = psi_out(x, spin, color);
              }
              dft(u.get(), v.get(), psi_out.lattice().size(0), sign);
              for (mdp_int i = 0; i < psi_out.lattice().size(0); i++)
              {
                x.set(i, x1, x2, x3);
                psi_out(x, spin, color) = u[i];
              }
            }
      }
  }

  /** @brief Default FT in space x-y-z
   *
   * Set ttime=true to FT in time too
   */
  void fermi_field_fft(fermi_field &psi_out,
                       fermi_field &psi_in,
                       int sign, bool ttime = false)
  {
    for (mdp_int t = 0; t < psi_in.lattice().size(0); t++)
      fermi_field_fft(t, psi_out, psi_in, sign);

    if (ttime)
      fermi_field_fft_t(psi_out, psi_out, sign);
  }

  void mdp_complex_field_fft(int t,
                             mdp_field<mdp_complex> &psi_out,
                             mdp_field<mdp_complex> &psi_in,
                             int sign)
  {

    if (psi_in.lattice().n_dimensions() != 4)
      error("fft3D requires TxXxXxX");

    mdp_int size = psi_in.lattice().size(1);
    if (psi_in.lattice().size(2) > size)
      size = psi_in.lattice().size(2);
    if (psi_in.lattice().size(3) > size)
      size = psi_in.lattice().size(3);

    std::unique_ptr<mdp_complex[]> v = std::make_unique<mdp_complex[]>(size);
    std::unique_ptr<mdp_complex[]> u = std::make_unique<mdp_complex[]>(size);

    mdp_site x(psi_in.lattice());

    forallsites(x)
    {
      if (x(0) == t)
        for (mdp_int k = 0; k < psi_in.size_per_site(); k++)
          psi_out(x, k) = psi_in(x, k);
    }

    for (mdp_int k = 0; k < psi_in.size_per_site(); k++)
    {
      for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
        for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
        {
          for (mdp_int i = 0; i < psi_out.lattice().size(1); i++)
          {
            x.set(t, i, x2, x3);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(1), sign);
          for (mdp_int i = 0; i < psi_out.lattice().size(1); i++)
          {
            x.set(t, i, x2, x3);
            psi_out(x, k) = u[i];
          }
        }
      for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
        {
          for (mdp_int i = 0; i < psi_out.lattice().size(2); i++)
          {
            x.set(t, x1, i, x3);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(2), sign);
          for (mdp_int i = 0; i < psi_out.lattice().size(2); i++)
          {
            x.set(t, x1, i, x3);
            psi_out(x, k) = u[i];
          }
        }

      for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
        {
          for (mdp_int i = 0; i < psi_out.lattice().size(3); i++)
          {
            x.set(t, x1, x2, i);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(3), sign);
          for (mdp_int i = 0; i < psi_out.lattice().size(3); i++)
          {
            x.set(t, x1, x2, i);
            psi_out(x, k) = u[i];
          }
        }
    }
  }

  void mdp_complex_field_fft_t(mdp_field<mdp_complex> &psi_out,
                               mdp_field<mdp_complex> &psi_in,
                               int sign)
  {

    if (psi_in.lattice().n_dimensions() != 4)
      error("fft3D requires TxXxXxX");

    mdp_int size = psi_in.lattice().size(0);
    if (psi_in.lattice().size(2) > size)
      size = psi_in.lattice().size(2);
    if (psi_in.lattice().size(3) > size)
      size = psi_in.lattice().size(3);

    std::unique_ptr<mdp_complex[]> v = std::make_unique<mdp_complex[]>(size);
    std::unique_ptr<mdp_complex[]> u = std::make_unique<mdp_complex[]>(size);

    mdp_site x(psi_in.lattice());

    forallsites(x)
    {
      for (mdp_int k = 0; k < psi_in.size_per_site(); k++)
        psi_out(x, k) = psi_in(x, k);
    }

    for (mdp_int k = 0; k < psi_in.size_per_site(); k++)
    {
      for (mdp_int x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (mdp_int x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          for (mdp_int x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (mdp_int i = 0; i < psi_out.lattice().size(0); i++)
            {
              x.set(i, x1, x2, x3);
              v[i] = psi_out(x, k);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(0), sign);
            for (mdp_int i = 0; i < psi_out.lattice().size(0); i++)
            {
              x.set(i, x1, x2, x3);
              psi_out(x, k) = u[i];
            }
          }
    }
  }

  /** @brief Default FT in space x-y-z
   *
   * Set ttime=true to FT in time too
   */
  void mdp_complex_field_fft(mdp_field<mdp_complex> &psi_out,
                             mdp_field<mdp_complex> &psi_in,
                             int sign, bool ttime = false)
  {
    for (mdp_int t = 0; t < psi_in.lattice().size(0); t++)
      mdp_complex_field_fft(t, psi_out, psi_in, sign);

    if (ttime)
      mdp_complex_field_fft_t(psi_out, psi_out, sign);
  }
} // namespace MDP

#endif /* FERMIQCD_FFTS_ */
