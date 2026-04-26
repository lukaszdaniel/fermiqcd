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
  /** @brief Computes the one-dimensional Discrete Fourier Transform (DFT)
   *
   * This function computes the DFT (or inverse DFT) of a complex input array `f`
   * and stores the result in `fft_f`.
   *
   * The transform is defined as:
   *   fft_f[i] = (1 / sqrt(n)) * sum_{j=0}^{n-1} f[j] * exp(2 * Pi * I * sign * i * j / n)
   *
   * where:
   *   - sign = -1 -- forward DFT
   *   - sign = +1 -- inverse DFT
   *
   * Parameters `offset` and `coeff` allow operating on strided data, enabling
   * this function to be used on subarrays or multi-dimensional data layouts.
   *
   * @param fft_f  Output array (size at least offset + coeff * (n-1))
   * @param f      Input array (size at least offset + coeff * (n-1))
   * @param n      Number of elements in the transform (DFT size)
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   * @param offset Starting index in the arrays
   * @param coeff  Stride between consecutive elements
   *
   * @note This is a naive O(n^2) implementation (not FFT).
   * @note Uses symmetric normalization (1 / sqrt(n)).
   */
  void dft(std::vector<mdp_complex> &fft_f, const std::vector<mdp_complex> &f, mdp_uint n, mdp_sint sign,
           mdp_uint offset = 0, mdp_uint coeff = 1)
  {
    if (sign == 0)
    {
      for (mdp_uint i = 0; i < n; i++)
        fft_f[offset + coeff * i] = f[offset + coeff * i];
      return;
    }

    mdp_complex phase = exp(mdp_complex(0.0, sign * 2.0 * Pi / n));
    mdp_complex w_i = 1.0f;

    for (mdp_uint i = 0; i < n; i++)
    {
      mdp_complex sum = 0.0f;
      mdp_complex w = 1.0f;

      for (mdp_uint j = 0; j < n; j++)
      {
        sum += f[offset + coeff * j] * w;
        w *= w_i;
      }

      fft_f[offset + coeff * i] = sum / std::sqrt((float)n);

      w_i *= phase;
    }
  }

  /**
   * @brief Computes the one-dimensional Fast Fourier Transform (FFT)
   *        using the iterative radix-2 Cooley–Tukey algorithm.
   *
   * This function computes the discrete Fourier transform (DFT) or its inverse
   * for an input array `f` of size N = 2^n and stores the result in `fft_f`.
   *
   * The transform is defined as:
   *   fft_f[i] = (1 / sqrt(N)) * sum_{j=0}^{N-1} f[j] * exp(2 * Pi * I * sign * i * j / N)
   *
   * where:
   *   - sign = -1 -- forward FFT
   *   - sign = +1 -- inverse FFT
   *   - sign =  0 -- no transform (copy input to output)
   *
   * The implementation uses:
   *   - bit-reversal permutation
   *   - iterative Cooley–Tukey radix-2 decomposition
   *
   * Parameters `offset` and `coeff` allow operating on strided data, enabling
   * this function to be used on subarrays or multi-dimensional data layouts.
   *
   * @param fft_f  Output array (size at least offset + coeff * (N-1))
   * @param f      Input array (size at least offset + coeff * (N-1))
   * @param n      Number of elements in the transform (must be a power of two)
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   * @param offset Starting index in the arrays
   * @param coeff  Stride between consecutive elements
   *
   * @note Time complexity: O(N log N)
   * @note Uses symmetric normalization (1 / sqrt(N))
   * @note Input size must be a power of two
   */
  void fft(std::vector<mdp_complex> &fft_f, const std::vector<mdp_complex> &f, mdp_uint n, mdp_sint sign,
           mdp_uint offset = 0, mdp_uint coeff = 1)
  {
    if (!n || (n & (n - 1)))
      error("fft requires n to be a power of two");

    if (sign == 0)
    {
      for (mdp_uint i = 0; i < n; i++)
        fft_f[offset + coeff * i] = f[offset + coeff * i];
      return;
    }

    // Copy data from f to fft_f
    for (mdp_uint i = 0; i < n; i++)
      fft_f[offset + coeff * i] = f[offset + coeff * i];

    // Bit reversal step (important for FFT)
    mdp_uint j = 0;
    for (mdp_uint i = 1; i < n; i++)
    {
      mdp_uint bit = n >> 1;
      for (; j >= bit; bit >>= 1)
        j -= bit;
      j += bit;

      if (i < j)
        std::swap(fft_f[offset + coeff * i], fft_f[offset + coeff * j]);
    }

    // Cooley-Tukey FFT
    for (mdp_uint m = 2; m <= n; m <<= 1)
    {
      mdp_uint half_m = m >> 1;
      mdp_complex omega_m = exp(mdp_complex(0.0, sign * 2.0 * Pi / m)); // Twiddle factor

      for (mdp_uint k = 0; k < n; k += m)
      {
        mdp_complex omega = 1.0;

        for (mdp_uint j = 0; j < half_m; j++)
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
    for (mdp_uint i = 0; i < n; i++)
    {
      fft_f[offset + coeff * i] /= std::sqrt(mdp_real(n));
    }
  }

  /**
   * @brief Hybrid DFT/FFT (automatic selection)
   *
   * Uses FFT when possible (power-of-two size and contiguous data),
   * otherwise falls back to DFT.
   */
  void ft_auto(std::vector<mdp_complex> &fft_f,
               const std::vector<mdp_complex> &f,
               mdp_uint n,
               mdp_sint sign,
               mdp_uint offset = 0,
               mdp_uint coeff = 1)
  {
    if (sign == 0)
    {
      for (mdp_uint i = 0; i < n; i++)
        fft_f[offset + coeff * i] = f[offset + coeff * i];
      return;
    }

    auto is_power_of_two = [](mdp_uint n)
    { return n && !(n & (n - 1)); };

    if (coeff == 1 && is_power_of_two(n))
    {
      fft(fft_f, f, n, sign, offset, coeff);
    }
    else
    {
      dft(fft_f, f, n, sign, offset, coeff);
    }
  }

  /**
   * @brief Computes a 3D Discrete Fourier Transform (DFT) on a fixed time slice
   *        of a 4D fermion field.
   *
   * This function computes the DFT independently along spatial dimensions
   * (x, y, z) for a given time slice `t` of the input field `psi_in`.
   * The result is stored in `psi_out`.
   *
   * The transform is applied separately for each spin and color component.
   *
   * where:
   *   - sign = -1 -- forward FFT
   *   - sign = +1 -- inverse FFT
   *   - sign =  0 -- no transform (copy input to output)
   *
   * The implementation uses separability of the Fourier transform:
   *   DFT_3D = DFT_x * DFT_y * DFT_z
   *
   * @param t        Time slice index
   * @param psi_out  Output fermion field
   * @param psi_in   Input fermion field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   *
   * @note Only the slice x(0) == t is transformed; other slices are untouched.
   * @note This is a naive DFT implementation with O(N^2) complexity per axis.
   * @note Uses 1D DFT applied successively along each spatial dimension.
   */
  void fermi_field_fft(mdp_uint t,
                       fermi_field &psi_out,
                       const fermi_field &psi_in,
                       mdp_sint sign)
  {
    const auto &lat = psi_in.lattice();

    if (lat.n_dimensions() != 4)
      error("fermi_field_fft requires 4D lattice of form TxXxXxX");

    const mdp_uint Nx = lat.size(1);
    const mdp_uint Ny = lat.size(2);
    const mdp_uint Nz = lat.size(3);

    const mdp_uint max_size = std::max({Nx, Ny, Nz});

    std::vector<mdp_complex> v(max_size);
    std::vector<mdp_complex> u(max_size);

    mdp_site x(lat);

    if (&psi_out != &psi_in)
    {
      forallsites(x)
      {
        if (x(0) == t)
          psi_out(x) = psi_in(x);
      }
    }

    for (mdp_suint spin = 0; spin < psi_out.nspin(); spin++)
      for (mdp_suint color = 0; color < psi_out.nc(); color++)
      {
        for (mdp_uint y = 0; y < Ny; y++)
          for (mdp_uint z = 0; z < Nz; z++)
          {
            for (mdp_uint i = 0; i < Nx; i++)
            {
              x.set(t, i, y, z);
              v[i] = psi_out(x, spin, color);
            }

            ft_auto(u, v, Nx, sign);

            for (mdp_uint i = 0; i < Nx; i++)
            {
              x.set(t, i, y, z);
              psi_out(x, spin, color) = u[i];
            }
          }

        for (mdp_uint x1 = 0; x1 < Nx; x1++)
          for (mdp_uint z = 0; z < Nz; z++)
          {
            for (mdp_uint i = 0; i < Ny; i++)
            {
              x.set(t, x1, i, z);
              v[i] = psi_out(x, spin, color);
            }

            ft_auto(u, v, Ny, sign);

            for (mdp_uint i = 0; i < Ny; i++)
            {
              x.set(t, x1, i, z);
              psi_out(x, spin, color) = u[i];
            }
          }

        for (mdp_uint x1 = 0; x1 < Nx; x1++)
          for (mdp_uint y = 0; y < Ny; y++)
          {
            for (mdp_uint i = 0; i < Nz; i++)
            {
              x.set(t, x1, y, i);
              v[i] = psi_out(x, spin, color);
            }

            ft_auto(u, v, Nz, sign);

            for (mdp_uint i = 0; i < Nz; i++)
            {
              x.set(t, x1, y, i);
              psi_out(x, spin, color) = u[i];
            }
          }
      }
  }

  /**
   * @brief Computes a 1D Discrete Fourier Transform (DFT) along the time dimension
   *        of a 4D fermion field.
   *
   * This function computes the DFT independently along the time direction (t)
   * for each spatial point (x, y, z) of the input field `psi_in`.
   * The result is stored in `psi_out`.
   *
   * The transform is applied separately for each spin and color component.
   *
   * where:
   *   - sign = -1 -- forward FFT
   *   - sign = +1 -- inverse FFT
   *   - sign =  0 -- no transform (copy input to output)
   *
   * @param psi_out  Output fermion field
   * @param psi_in   Input fermion field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   *
   * @note This is a naive DFT implementation with O(N^2) complexity per axis.
   * @note Uses 1D DFT applied successively along each spatial dimension.
   */
  void fermi_field_fft_t(fermi_field &psi_out,
                         const fermi_field &psi_in,
                         mdp_sint sign)
  {
    const auto &lat = psi_in.lattice();

    if (lat.n_dimensions() != 4)
      error("fermi_field_fft_t requires 4D lattice of form TxXxXxX");

    const mdp_uint Nt = lat.size(0);
    const mdp_uint Nx = lat.size(1);
    const mdp_uint Ny = lat.size(2);
    const mdp_uint Nz = lat.size(3);

    const mdp_uint max_size = std::max({Nt, Ny, Nz});

    std::vector<mdp_complex> v(max_size);
    std::vector<mdp_complex> u(max_size);

    mdp_site x(lat);

    if (&psi_out != &psi_in)
    {
      forallsites(x)
      {
        psi_out(x) = psi_in(x);
      }
    }

    for (mdp_suint spin = 0; spin < psi_out.nspin(); spin++)
      for (mdp_suint color = 0; color < psi_out.nc(); color++)
      {
        for (mdp_uint x1 = 0; x1 < Nx; x1++)
          for (mdp_uint y = 0; y < Ny; y++)
            for (mdp_uint z = 0; z < Nz; z++)
            {
              for (mdp_uint t = 0; t < Nt; t++)
              {
                x.set(t, x1, y, z);
                v[t] = psi_out(x, spin, color);
              }

              ft_auto(u, v, Nt, sign);

              for (mdp_uint t = 0; t < Nt; t++)
              {
                x.set(t, x1, y, z);
                psi_out(x, spin, color) = u[t];
              }
            }
      }
  }

  /**
   * @brief Computes a 3+1 Discrete Fourier Transform (DFT) of a 4D fermion field.
   *
   * By default, applies 3D DFT in spatial dimensions (x, y, z)
   * independently for each time slice.
   *
   * If ttime = true, an additional DFT is applied in the time dimension,
   * resulting in a full 4D Fourier transform.
   *
   * @param psi_out  Output fermion field
   * @param psi_in   Input fermion field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   * @param ttime   If true, also perform DFT in time direction
   */
  void fermi_field_fft(fermi_field &psi_out,
                       const fermi_field &psi_in,
                       mdp_sint sign, bool ttime = false)
  {
    const mdp_uint Nt = psi_in.lattice().size(0);

    for (mdp_uint t = 0; t < Nt; t++)
      fermi_field_fft(t, psi_out, psi_in, sign);

    if (ttime)
      fermi_field_fft_t(psi_out, psi_out, sign);
  }

  /**
   * @brief Computes a 3D Discrete Fourier Transform (DFT) on a fixed time slice
   *        of a 4D complex field.
   *
   * This function computes the DFT independently along spatial dimensions
   * (x, y, z) for a given time slice `t` of the input field `psi_in`.
   * The result is stored in `psi_out`.
   *
   * The transform is applied separately for each internal field component.
   *
   * where:
   *   - sign = -1 -- forward FFT
   *   - sign = +1 -- inverse FFT
   *   - sign =  0 -- no transform (copy input to output)
   *
   * The implementation uses separability of the Fourier transform:
   *   DFT_3D = DFT_x * DFT_y * DFT_z
   *
   * @param t        Time slice index
   * @param psi_out  Output complex field
   * @param psi_in   Input complex field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   *
   * @note Only the slice x(0) == t is transformed; other slices are untouched.
   * @note This is a naive DFT implementation with O(N^2) complexity per axis.
   * @note Uses 1D DFT applied successively along each spatial dimension.
   */
  void mdp_complex_field_fft(mdp_uint t,
                             mdp_complex_field &psi_out,
                             const mdp_complex_field &psi_in,
                             mdp_sint sign)
  {
    const auto &lat = psi_in.lattice();

    if (lat.n_dimensions() != 4)
      error("mdp_complex_field_fft requires 4D lattice of form TxXxXxX");

    const mdp_uint Nx = lat.size(1);
    const mdp_uint Ny = lat.size(2);
    const mdp_uint Nz = lat.size(3);

    const mdp_uint max_size = std::max({Nx, Ny, Nz});

    std::vector<mdp_complex> v(max_size);
    std::vector<mdp_complex> u(max_size);

    mdp_site x(lat);

    if (&psi_out != &psi_in)
    {
      forallsites(x)
      {
        if (x(0) == t)
          for (mdp_uint k = 0; k < psi_in.size_per_site(); k++)
            psi_out(x, k) = psi_in(x, k);
      }
    }

    for (mdp_uint k = 0; k < psi_in.size_per_site(); k++)
    {
      for (mdp_uint y = 0; y < Ny; y++)
        for (mdp_uint z = 0; z < Nz; z++)
        {
          for (mdp_uint i = 0; i < Nx; i++)
          {
            x.set(t, i, y, z);
            v[i] = psi_out(x, k);
          }

          ft_auto(u, v, Nx, sign);

          for (mdp_uint i = 0; i < Nx; i++)
          {
            x.set(t, i, y, z);
            psi_out(x, k) = u[i];
          }
        }

      for (mdp_uint x1 = 0; x1 < Nx; x1++)
        for (mdp_uint z = 0; z < Nz; z++)
        {
          for (mdp_uint i = 0; i < Ny; i++)
          {
            x.set(t, x1, i, z);
            v[i] = psi_out(x, k);
          }

          ft_auto(u, v, Ny, sign);

          for (mdp_uint i = 0; i < Ny; i++)
          {
            x.set(t, x1, i, z);
            psi_out(x, k) = u[i];
          }
        }

      for (mdp_uint x1 = 0; x1 < Nx; x1++)
        for (mdp_uint y = 0; y < Ny; y++)
        {
          for (mdp_uint i = 0; i < Nz; i++)
          {
            x.set(t, x1, y, i);
            v[i] = psi_out(x, k);
          }

          ft_auto(u, v, Nz, sign);

          for (mdp_uint i = 0; i < Nz; i++)
          {
            x.set(t, x1, y, i);
            psi_out(x, k) = u[i];
          }
        }
    }
  }

  /**
   * @brief Computes a 1D Discrete Fourier Transform (DFT) along the time dimension
   *        of a 4D complex field.
   *
   * This function computes the DFT independently along the time direction (t)
   * for each spatial point (x, y, z) of the input field `psi_in`.
   * The result is stored in `psi_out`.
   *
   * The transform is applied separately for each internal field component.
   *
   * where:
   *   - sign = -1 -- forward FFT
   *   - sign = +1 -- inverse FFT
   *   - sign =  0 -- no transform (copy input to output)
   *
   * @param psi_out  Output fermion field
   * @param psi_in   Input fermion field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   *
   * @note This is a naive DFT implementation with O(N^2) complexity per axis.
   * @note Uses 1D DFT applied successively along each spatial dimension.
   */
  void mdp_complex_field_fft_t(mdp_complex_field &psi_out,
                               const mdp_complex_field &psi_in,
                               mdp_sint sign)
  {
    const auto &lat = psi_in.lattice();

    if (lat.n_dimensions() != 4)
      error("mdp_complex_field_fft_t requires 4D lattice of form TxXxXxX");

    const mdp_uint Nt = lat.size(0);
    const mdp_uint Nx = lat.size(1);
    const mdp_uint Ny = lat.size(2);
    const mdp_uint Nz = lat.size(3);

    const mdp_uint max_size = std::max({Nt, Ny, Nz});

    std::vector<mdp_complex> v(max_size);
    std::vector<mdp_complex> u(max_size);

    mdp_site x(lat);

    if (&psi_out != &psi_in)
    {
      forallsites(x)
      {
        for (mdp_uint k = 0; k < psi_in.size_per_site(); k++)
          psi_out(x, k) = psi_in(x, k);
      }
    }

    for (mdp_uint k = 0; k < psi_in.size_per_site(); k++)
    {
      for (mdp_uint x1 = 0; x1 < Nx; x1++)
        for (mdp_uint y = 0; y < Ny; y++)
          for (mdp_uint z = 0; z < Nz; z++)
          {
            for (mdp_uint t = 0; t < Nt; t++)
            {
              x.set(t, x1, y, z);
              v[t] = psi_out(x, k);
            }

            ft_auto(u, v, Nt, sign);

            for (mdp_uint t = 0; t < Nt; t++)
            {
              x.set(t, x1, y, z);
              psi_out(x, k) = u[t];
            }
          }
    }
  }

  /**
   * @brief Computes Fourier transform of a 4D complex field.
   *
   * By default, applies 3D DFT in spatial dimensions (x, y, z)
   * independently for each time slice.
   *
   * If ttime = true, an additional DFT is applied in the time dimension,
   * resulting in a full 4D Fourier transform.
   *
   * @param psi_out  Output complex field
   * @param psi_in   Input complex field
   * @param sign   Direction of transform: -1 (forward), +1 (inverse), 0 (copy only)
   * @param ttime   If true, also perform DFT in time direction
   */
  void mdp_complex_field_fft(mdp_complex_field &psi_out,
                             const mdp_complex_field &psi_in,
                             mdp_sint sign, bool ttime = false)
  {
    const mdp_uint Nt = psi_in.lattice().size(0);

    for (mdp_uint t = 0; t < Nt; t++)
      mdp_complex_field_fft(t, psi_out, psi_in, sign);

    if (ttime)
      mdp_complex_field_fft_t(psi_out, psi_out, sign);
  }
} // namespace MDP

#endif /* FERMIQCD_FFTS_ */
