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
    mdp_complex phase = exp(2.0 * Pi * I * sign / n);
    for (int i = 0; i < n; i++)
    {
      fft_f[offset + coeff * i] = 0;
      for (int j = 0; j < n; j++)
        fft_f[offset + coeff * i] += f[offset + coeff * j] * pow(phase, i * j);
      fft_f[offset + coeff * i] /= std::sqrt(n);
    }
  }

#if 0
  // NOT NOT UNCOMMENT THIS, WORK IN PROGRESS!!!
  void fft(mdp_complex *fft_f, mdp_complex *f, mdp_int n, double sign,
           mdp_int offset = 0, mdp_int coeff = 1)
  {
    mdp_int a, b, h, pow2b, N = i2pow(n);
    mdp_complex alpha, omega, F[N];
    if (sign != 0)
      for (h = 0; h < N; h++)
      {
        alpha = 2.0 * sign * Pi * I / N * h;
        for (a = 0; a < N; a++)
          F[a] = f[offset + coeff * a];
        for (b = n - 1; b >= 0; b--)
        {
          pow2b = i2pow(b);
          omega = exp(alpha * pow2b);
          for (a = 0; a < pow2b; a++)
            F[a] += omega * F[a + pow2b];
        };
        fft_f[offset + coeff * h] = F[0] / sqrt(N);
      }
    else
    {
      for (a = 0; a < N; a++)
        fft_f[offset + coeff * a] = f[offset + coeff * a];
    }
  }
#endif

  void fermi_field_fft(int t,
                       fermi_field &psi_out,
                       fermi_field &psi_in,
                       int sign)
  {

    if (psi_in.lattice().n_dimensions() != 4)
      error("fft3D requires TxXxXxX");

    int i, x1, x2, x3, spin, color;
    int size = psi_in.lattice().size(1);
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

    for (spin = 0; spin < psi_out.nspin; spin++)
      for (color = 0; color < psi_out.nc; color++)
      {
        for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (i = 0; i < psi_out.lattice().size(1); i++)
            {
              x.set(t, i, x2, x3);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(1), sign);
            for (i = 0; i < psi_out.lattice().size(1); i++)
            {
              x.set(t, i, x2, x3);
              psi_out(x, spin, color) = u[i];
            }
          }
        for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (i = 0; i < psi_out.lattice().size(2); i++)
            {
              x.set(t, x1, i, x3);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(2), sign);
            for (i = 0; i < psi_out.lattice().size(2); i++)
            {
              x.set(t, x1, i, x3);
              psi_out(x, spin, color) = u[i];
            }
          }

        for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          {
            for (i = 0; i < psi_out.lattice().size(3); i++)
            {
              x.set(t, x1, x2, i);
              v[i] = psi_out(x, spin, color);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(3), sign);
            for (i = 0; i < psi_out.lattice().size(3); i++)
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

    int i, x1, x2, x3, spin, color;
    int size = psi_in.lattice().size(0);
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

    for (spin = 0; spin < psi_out.nspin; spin++)
      for (color = 0; color < psi_out.nc; color++)
      {
        for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
          for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
            for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
            {
              for (i = 0; i < psi_out.lattice().size(0); i++)
              {
                x.set(i, x1, x2, x3);
                v[i] = psi_out(x, spin, color);
              }
              dft(u.get(), v.get(), psi_out.lattice().size(0), sign);
              for (i = 0; i < psi_out.lattice().size(0); i++)
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
    for (int t = 0; t < psi_in.lattice().size(0); t++)
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

    int i, x1, x2, x3;
    int size = psi_in.lattice().size(1);
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
        for (int k = 0; k < psi_in.size_per_site(); k++)
          psi_out(x, k) = psi_in(x, k);
    }

    for (int k = 0; k < psi_in.size_per_site(); k++)
    {
      for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
        for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
        {
          for (i = 0; i < psi_out.lattice().size(1); i++)
          {
            x.set(t, i, x2, x3);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(1), sign);
          for (i = 0; i < psi_out.lattice().size(1); i++)
          {
            x.set(t, i, x2, x3);
            psi_out(x, k) = u[i];
          }
        }
      for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
        {
          for (i = 0; i < psi_out.lattice().size(2); i++)
          {
            x.set(t, x1, i, x3);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(2), sign);
          for (i = 0; i < psi_out.lattice().size(2); i++)
          {
            x.set(t, x1, i, x3);
            psi_out(x, k) = u[i];
          }
        }

      for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
        {
          for (i = 0; i < psi_out.lattice().size(3); i++)
          {
            x.set(t, x1, x2, i);
            v[i] = psi_out(x, k);
          }
          dft(u.get(), v.get(), psi_out.lattice().size(3), sign);
          for (i = 0; i < psi_out.lattice().size(3); i++)
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

    int i, x1, x2, x3, k;
    int size = psi_in.lattice().size(0);
    if (psi_in.lattice().size(2) > size)
      size = psi_in.lattice().size(2);
    if (psi_in.lattice().size(3) > size)
      size = psi_in.lattice().size(3);

    std::unique_ptr<mdp_complex[]> v = std::make_unique<mdp_complex[]>(size);
    std::unique_ptr<mdp_complex[]> u = std::make_unique<mdp_complex[]>(size);

    mdp_site x(psi_in.lattice());

    forallsites(x)
    {
      for (k = 0; k < psi_in.size_per_site(); k++)
        psi_out(x, k) = psi_in(x, k);
    }

    for (k = 0; k < psi_in.size_per_site(); k++)
    {
      for (x1 = 0; x1 < psi_out.lattice().size(1); x1++)
        for (x2 = 0; x2 < psi_out.lattice().size(2); x2++)
          for (x3 = 0; x3 < psi_out.lattice().size(3); x3++)
          {
            for (i = 0; i < psi_out.lattice().size(0); i++)
            {
              x.set(i, x1, x2, x3);
              v[i] = psi_out(x, k);
            }
            dft(u.get(), v.get(), psi_out.lattice().size(0), sign);
            for (i = 0; i < psi_out.lattice().size(0); i++)
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
    for (int t = 0; t < psi_in.lattice().size(0); t++)
      mdp_complex_field_fft(t, psi_out, psi_in, sign);

    if (ttime)
      mdp_complex_field_fft_t(psi_out, psi_out, sign);
  }
} // namespace MDP

#endif /* FERMIQCD_FFTS_ */
