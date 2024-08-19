/////////////////////////////////////////////////////////////////
/// @file mdp_complex_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_complex_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_COMPLEX_FIELD_
#define MDP_COMPLEX_FIELD_

#include <cstdio>
#include "mdp_global_vars.h"
#include "mdp_field.h"

namespace MDP
{
  bool mdp_write_double_as_float(FILE *fp,
                                 void *data,
                                 mdp_int psize,
                                 mdp_int header_size,
                                 mdp_int position,
                                 const mdp_lattice &lattice)
  {
    double *p = (double *)data;
    float *q = (float *)malloc(psize / 2);
    for (mdp_uint i = 0; i < psize / sizeof(double); i++)
      q[i] = p[i];
    if (fseek(fp, position * psize / 2 + header_size, SEEK_SET) ||
        fwrite(q, psize / 2, 1, fp) != 1)
      return false;
    free(q);
    return true;
  }

  bool mdp_read_double_as_float(FILE *fp,
                                void *data,
                                mdp_int psize,
                                mdp_int header_size,
                                mdp_int position,
                                const mdp_lattice &lattice)
  {
    double *p = (double *)data;
    float *q = (float *)malloc(psize / 2);
    if (fseek(fp, position * psize / 2 + header_size, SEEK_SET) ||
        fread(q, psize / 2, 1, fp) != 1)
      return false;
    for (mdp_uint i = 0; i < psize / sizeof(double); i++)
      p[i] = q[i];
    free(q);
    return true;
  }

  bool mdp_write_float_as_double(FILE *fp,
                                 void *data,
                                 mdp_int psize,
                                 mdp_int header_size,
                                 mdp_int position,
                                 const mdp_lattice &lattice)
  {
    float *p = (float *)data;
    double *q = (double *)malloc(psize * 2);
    for (mdp_uint i = 0; i < psize / sizeof(float); i++)
      q[i] = p[i];
    if (fseek(fp, position * psize * 2 + header_size, SEEK_SET) ||
        fwrite(q, psize * 2, 1, fp) != 1)
      return false;
    free(q);
    return true;
  }

  bool mdp_read_float_as_double(FILE *fp,
                                void *data,
                                mdp_int psize,
                                mdp_int header_size,
                                mdp_int position,
                                const mdp_lattice &lattice)
  {
    float *p = (float *)data;
    double *q = (double *)malloc(psize * 2);
    if (fseek(fp, position * psize * 2 + header_size, SEEK_SET) ||
        fread(q, psize * 2, 1, fp) != 1)
      return false;
    for (mdp_uint i = 0; i < psize / sizeof(float); i++)
      p[i] = q[i];
    free(q);
    return true;
  }

  /// @brief field of complex numbers or vectors of complex numbers
  ///
  /// Example:
  /// @verbatim
  ///    int box[]={10,10,10};
  ///    mdp_lattice lattice(3,box);
  ///    mdp_complex_field psi(lattice,10);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      for(int i=0; i<10; i++)
  ///         psi(x,i)=0.0+0.0*I;
  /// @endverbatim
  class mdp_complex_field : public mdp_field<mdp_complex>
  {
  public:
    mdp_complex_field() : mdp_field<mdp_complex>()
    {
    }

    mdp_complex_field(mdp_lattice &lattice, int n = 1) : mdp_field<mdp_complex>(lattice, n)
    {
    }

    mdp_complex_field(const mdp_complex_field &other) : mdp_field<mdp_complex>(other)
    {
    }

    void operator=(const mdp_complex_field &psi)
    {
      if (&lattice() != &psi.lattice() ||
          m_size != psi.m_size ||
          m_field_components != psi.m_field_components)
        error("mdp_field: operator=() incompatible fields");

      mdp_uint i = 0;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
      _sse_double *r = (_sse_double *)m_data.get();
      _sse_double *s = (_sse_double *)psi.m_data.get();

      for (; i < m_size - 7; i += 8, r += 8, s += 8)
      {
        _sse_double_prefetch_16(s + 8);
        _sse_double_copy_16(r, s);
      }
#endif

      for (; i < m_size; i++)
        m_data[i] = psi.m_data[i];
    }

    void operator*=(const mdp_complex alpha)
    {
      mdp_int i_min = physical_local_start(EVENODD);
      mdp_int i_max = physical_local_stop(EVENODD);
      mdp_int i = i_min;
      for (; i < i_max; i++)
        m_data[i] *= alpha;
    }
    void operator/=(const mdp_complex alpha)
    {
      (*this) *= (1.0 / alpha);
    }

    void operator*=(const mdp_real alpha)
    {
      mdp_int i_min = physical_local_start(EVENODD);
      mdp_int i_max = physical_local_stop(EVENODD);
      mdp_int i = i_min;

      if (alpha == 1)
        return;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
      static _sse_double c ALIGN16;
      _sse_double *r = (_sse_double *)physical_address(i_min);

      _sse_check_alignment(&c, 0xf);

      c.c1 = c.c2 = alpha;
      for (; i < i_max - 7; i += 8, r += 8)
      {
        _sse_double_prefetch_16(r + 8);
        _sse_double_multiply_16(r, c, r);
      }
#endif

      for (; i < i_max; i++)
        m_data[i] *= alpha;
    }

    void operator/=(const mdp_real alpha)
    {
      (*this) *= (1.0 / alpha);
    }

    void operator+=(mdp_complex_field &psi)
    {
      mdp_int i_min = psi.physical_local_start(EVENODD);
      mdp_int i_max = psi.physical_local_stop(EVENODD);
      mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
      _sse_double *r = (_sse_double *)physical_address(i_min);
      _sse_double *s = (_sse_double *)psi.physical_address(i_min);

      for (; i < i_max - 7; i += 8, r += 8, s += 8)
      {
        _sse_double_prefetch_16(s + 8);
        _sse_double_add_16(r, s);
      }
#endif

      for (; i < i_max; i++)
        m_data[i] += psi.m_data[i];
    }

    void operator-=(mdp_complex_field &psi)
    {
      mdp_int i_min = psi.physical_local_start(EVENODD);
      mdp_int i_max = psi.physical_local_stop(EVENODD);
      mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
      _sse_double *r = (_sse_double *)physical_address(i_min);
      _sse_double *s = (_sse_double *)psi.physical_address(i_min);

      for (; i < i_max - 7; i += 8, r += 8, s += 8)
      {
        _sse_double_prefetch_16(s + 8);
        _sse_double_sub_16(r, s);
      }
#endif

      for (; i < i_max; i++)
        m_data[i] -= psi.m_data[i];
    }

    bool save_as_float(std::string filename,
                       int processIO = 0,
                       mdp_int max_buffer_size = 1024,
                       bool load_header = true,
                       mdp_int skip_bytes = 0)
    {
#ifdef USE_DOUBLE_PRECISION
      m_header.bytes_per_site /= 2;
      save(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_write_double_as_float);
      m_header.bytes_per_site *= 2;
#else
      save(filename, processIO, max_buffer_size, load_header, skip_bytes, 0);
#endif
      return true;
    }

    bool load_as_float(std::string filename,
                       int processIO = 0,
                       mdp_int max_buffer_size = 1024,
                       bool load_header = true,
                       mdp_int skip_bytes = 0)
    {

#ifdef USE_DOUBLE_PRECISION
      m_header.bytes_per_site /= 2;
      load(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_read_double_as_float, true);
      m_header.bytes_per_site *= 2;
#else
      load(filename, processIO, max_buffer_size, load_header, skip_bytes, 0, true);
#endif
      return true;
    }

    bool load_as_double(std::string filename,
                        int processIO = 0,
                        mdp_int max_buffer_size = 1024,
                        bool load_header = true,
                        mdp_int skip_bytes = 0)
    {
#ifndef USE_DOUBLE_PRECISION
      m_header.bytes_per_site *= 2;
      load(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_read_float_as_double, true);
      m_header.bytes_per_site /= 2;
#else
      load(filename, processIO, max_buffer_size, load_header, skip_bytes, 0, true);
#endif
      return true;
    }

    bool save_as_double(std::string filename,
                        int processIO = 0,
                        mdp_int max_buffer_size = 1024,
                        bool load_header = true,
                        mdp_int skip_bytes = 0)
    {
#ifndef USE_DOUBLE_PRECISION
      m_header.bytes_per_site *= 2;
      save(filename, processIO, max_buffer_size, load_header, skip_bytes,
           mdp_write_float_as_double);
      m_header.bytes_per_site /= 2;
#else
      save(filename, processIO, max_buffer_size, load_header, skip_bytes, 0);
#endif
      return true;
    }
  };

  mdp_real norm_square(mdp_complex_field &psi,
                       int parity = EVENODD)
  {
    double n2 = 0;
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);
    mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
    static _sse_double c ALIGN16;
    _sse_double *r = (_sse_double *)psi.physical_address(i_min);

    _sse_check_alignment(&c, 0xf);

    c.c1 = c.c2 = 0;
    for (; i < i_max - 7; i += 8, r += 8)
    {
      _sse_double_prefetch_16(r + 8);
      _sse_double_add_norm_square_16(r, c);
    }
    n2 += c.c1 + c.c2;
#endif

    for (; i < i_max; i++)
      n2 += abs2(psi[i]);
    mdp.add(n2);
    return n2;
  }

  mdp_complex scalar_product(mdp_complex_field &psi,
                             mdp_complex_field &chi,
                             int parity = EVENODD)
  {
    mdp_complex n2 = 0;
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);
    mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
    static _sse_double c, d ALIGN16;
    _sse_double *r = (_sse_double *)psi.physical_address(i_min);
    _sse_double *s = (_sse_double *)chi.physical_address(i_min);

    _sse_check_alignment(&c, 0xf);

    c.c1 = c.c2 = 0;
    d.c1 = d.c2 = 0;
    for (; i < i_max - 7; i += 8, r += 8, s += 8)
    {
      _sse_double_prefetch_16(r + 8);
      _sse_double_prefetch_16(s + 8);
      _sse_double_add_real_scalar_product_16(r, s, c);
      _sse_double_add_imag_scalar_product_16(r, s, d);
    }
    n2 += mdp_complex(c.c1 + c.c2, d.c2 - d.c1);
#endif

    for (; i < i_max; i++)
      n2 += conj(psi[i]) * chi[i];
    mdp.add(n2);

    return n2;
  }

  mdp_real real_scalar_product(mdp_complex_field &psi,
                               mdp_complex_field &chi,
                               int parity = EVENODD)
  {

    double n2 = 0;
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);
    mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
    static _sse_double c ALIGN16;
    _sse_double *r = (_sse_double *)psi.physical_address(i_min);
    _sse_double *s = (_sse_double *)chi.physical_address(i_min);

    _sse_check_alignment(&c, 0xf);

    c.c1 = c.c2 = 0;
    for (; i < i_max - 7; i += 8, r += 8, s += 8)
    {
      _sse_double_prefetch_16(r + 8);
      _sse_double_prefetch_16(s + 8);
      _sse_double_add_real_scalar_product_16(r, s, c);
    }
    n2 += c.c1 + c.c2;
#endif

    for (; i < i_max; i++)
      n2 +=
          real(chi[i]) * real(psi[i]) +
          imag(chi[i]) * imag(psi[i]);

    mdp.add(n2);
    return n2;
  }

  mdp_real imag_scalar_product(mdp_complex_field &psi,
                               mdp_complex_field &chi,
                               int parity = EVENODD)
  {
    double n2 = 0;
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);
    mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
    static _sse_double c ALIGN16;
    _sse_double *r = (_sse_double *)psi.physical_address(i_min);
    _sse_double *s = (_sse_double *)chi.physical_address(i_min);

    _sse_check_alignment(&c, 0xf);

    c.c1 = c.c2 = 0;
    for (; i < i_max - 7; i += 8, r += 8, s += 8)
    {
      _sse_double_prefetch_16(r + 8);
      _sse_double_prefetch_16(s + 8);
      _sse_double_add_imag_scalar_product_16(r, s, c);
    }
    n2 += c.c2 - c.c1;
#endif

    for (; i < i_max; i++)
      n2 +=
          real(psi[i]) * imag(chi[i]) +
          imag(psi[i]) * real(chi[i]);
    mdp.add(n2);
    return n2;
  }

  void mdp_add_scaled_field(mdp_complex_field &psi,
                            mdp_real alpha,
                            mdp_complex_field &chi,
                            int parity = EVENODD)
  {
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);
    mdp_int i = i_min;

#if defined(SSE2) && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
    static _sse_double c ALIGN16;
    _sse_double *r = (_sse_double *)psi.physical_address(i_min);
    _sse_double *s = (_sse_double *)chi.physical_address(i_min);

    _sse_check_alignment(&c, 0xf);

    c.c1 = c.c2 = alpha;
    for (i = 0; i < i_max - 7; i += 8, r += 8, s += 8)
    {
      _sse_double_prefetch_16(r + 8);
      _sse_double_prefetch_16(s + 8);
      _sse_double_add_multiply_16(r, c, s);
    }

#endif

    for (; i < i_max; i++)
      psi[i] += alpha * chi[i];
  }

  void mdp_add_scaled_field(mdp_complex_field &psi,
                            mdp_complex alpha,
                            mdp_complex_field &chi,
                            int parity = EVENODD)
  {
    mdp_int i_min = psi.physical_local_start(parity);
    mdp_int i_max = psi.physical_local_stop(parity);

    //    this needs optimization.
    for (mdp_int i = i_min; i < i_max; i++)
      psi[i] += alpha * chi[i];
  }

  mdp_complex operator*(mdp_complex_field &psi,
                        mdp_complex_field &chi)
  {
    return scalar_product(psi, chi);
  }

  mdp_real relative_residue(mdp_complex_field &p,
                            mdp_complex_field &q,
                            int parity = EVENODD)
  {
    double residue = 0, num = 0, den = 0;
    mdp_int i_min = p.physical_local_start(parity);
    mdp_int i_max = q.physical_local_stop(parity);

    //    this needs optimization.
    for (mdp_int i = i_min; i < i_max;)
    {
      num += abs2(p[i]);
      den += abs2(q[i]);
      if (++i % p.size_per_site() == 0)
      {
        residue += (den == 0) ? 1.0 : (num / den);
        num = den = 0;
      }
    }
    mdp.add(residue);
    return std::sqrt(residue / p.lattice().global_volume());
  }
} // namespace MDP

#endif /* MDP_COMPLEX_FIELD_ */
