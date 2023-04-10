/////////////////////////////////////////////////////////////////
/// @file mdp_matrix.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_matrix
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_MATRIX_
#define MDP_MATRIX_

#include <iostream>
#include <memory>

namespace MDP
{
  /// @brief matrices of complex numbers
  ///
  /// Example:
  /// @verbatim
  ///    mdp_matrix A,B;
  ///    A.dimension(3,3);
  ///    A(0,0)=A(1,1)=A(2,2)=A(1,2)=1.0+I/2;
  ///    B=A+inv(A)+exp(A+5);
  /// @endverbatim
  class mdp_matrix
  {
  private:
    bool m_shared;
    mdp_uint m_rows;
    mdp_uint m_cols;
    std::unique_ptr<mdp_complex[]> m_data; // data

    /** @brief Allocate required memory
     */
    void allocate()
    {
      if (!m_shared)
      {
        m_data = std::make_unique<mdp_complex[]>(size());
        if (size() && !m_data)
          error("mdp_matrix::allocate()\nOut of memory");
        // memset(m_data.get(), 0, size() * sizeof(mdp_complex));
      }
    }

    void reallocate()
    {
      deallocate();
      allocate();
    }

    void deallocate()
    {
      if (m_shared)
      {
        m_data.release();
      }
      else
      {
        m_data = nullptr;
      }
    }

  public:
    void dimension(mdp_uint r, mdp_uint c)
    {
      m_shared = false;
      m_rows = r;
      m_cols = c;
      reallocate();
    }

    mdp_matrix() : m_shared(false), m_rows(0), m_cols(0), m_data(nullptr)
    {
      allocate();
    }

    mdp_matrix(const mdp_matrix &a) : m_shared(false), m_rows(a.m_rows), m_cols(a.m_cols), m_data(nullptr)
    {
      allocate();
      for (mdp_uint i = 0; i < size(); i++)
        m_data[i] = a[i];
    }

    mdp_matrix(mdp_uint r, mdp_uint c) : m_shared(false), m_rows(r), m_cols(c), m_data(nullptr)
    {
      allocate();
    }

    mdp_matrix(mdp_complex *z, mdp_uint r, mdp_uint c) : m_shared(true), m_rows(r), m_cols(c), m_data(z)
    {
#if defined(MATRIX_SSE2) && defined(USE_DOUBLE_PRECISION)
      _sse_check_alignment((void *)m_data, 0xf);
#endif
    }

    virtual ~mdp_matrix()
    {
      deallocate();
    }

    const mdp_matrix &operator=(const mdp_matrix &x)
    {
      if (rows() != x.rows() || cols() != x.cols())
      {
        m_rows = x.rows();
        m_cols = x.cols();
        reallocate();
      }

      for (mdp_uint i = 0; i < size(); i++)
        m_data[i] = x[i];
      return *this;
    }

    mdp_complex &operator()(mdp_uint i, mdp_uint j)
    {
      return m_data[i * cols() + j];
    }

    const mdp_complex &operator()(mdp_uint i, mdp_uint j) const
    {
      return m_data[i * cols() + j];
    }

    mdp_matrix operator()(mdp_uint i)
    {
      return mdp_matrix(m_data.get() + i * cols(), cols(), 1);
    }

    /** @brief Begining of this matrix
     *
     * @return Pointer to the begining of this matrix
     */
    mdp_complex *address()
    {
      return m_data.get();
    }

    mdp_complex &operator[](mdp_uint i)
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_array::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    const mdp_complex &operator[](mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_array::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    /** @brief Size of i-th dimension
     *
     * @return Size of i-th dimension
     */
    mdp_uint rows() const
    {
      return m_rows;
    }

    /** @brief Size of i-th dimension
     *
     * @return Size of i-th dimension
     */
    mdp_uint cols() const
    {
      return m_cols;
    }

    /** @brief Total size of this array
     *
     * @return Total number of elements
     *
     * @note This is a product of all dimensions
     */
    mdp_uint size() const
    {
      return m_rows * m_cols;
    }

    mdp_matrix operator+=(const mdp_matrix &a)
    {
      (*this) = (*this) + a;
      return *this;
    }

    mdp_matrix operator-=(const mdp_matrix &a)
    {
      (*this) = (*this) - a;
      return *this;
    }

    mdp_matrix operator*=(const mdp_matrix &a)
    {
      (*this) = (*this) * a;
      return *this;
    }

    mdp_matrix operator/=(const mdp_matrix &a)
    {
      (*this) = (*this) / a;
      return *this;
    }

    mdp_matrix operator+=(mdp_complex a)
    {
      (*this) = (*this) + a;
      return *this;
    }

    mdp_matrix operator-=(mdp_complex a)
    {
      (*this) = (*this) - a;
      return *this;
    }

    mdp_matrix operator*=(mdp_complex a)
    {
      (*this) = (*this) * a;
      return *this;
    }

    mdp_matrix operator/=(mdp_complex a)
    {
      (*this) = (*this) / a;
      return *this;
    }

    mdp_matrix operator+=(mdp_real a)
    {
      (*this) = (*this) + a;
      return *this;
    }

    mdp_matrix operator-=(mdp_real a)
    {
      (*this) = (*this) - a;
      return *this;
    }

    mdp_matrix operator*=(mdp_real a)
    {
      (*this) = (*this) * a;
      return *this;
    }

    mdp_matrix operator/=(mdp_real a)
    {
      (*this) = (*this) / a;
      return *this;
    }

    void operator=(mdp_complex a)
    {
      for (mdp_uint i = 0; i < rows(); i++)
        for (mdp_uint j = 0; j < cols(); j++)
          (*this)(i, j) = (i == j) ? a : 0;
    }

    void operator=(mdp_real a)
    {
      for (mdp_uint i = 0; i < rows(); i++)
        for (mdp_uint j = 0; j < cols(); j++)
          (*this)(i, j) = (i == j) ? mdp_complex(a, 0) : 0;
    }

    friend mdp_matrix operator+(const mdp_matrix &a);
    friend mdp_matrix operator-(const mdp_matrix &a);

    friend mdp_matrix operator+(const mdp_matrix &a, const mdp_matrix &b);
    friend mdp_matrix operator-(const mdp_matrix &a, const mdp_matrix &b); // sse2
    friend mdp_matrix operator*(const mdp_matrix &a, const mdp_matrix &b); // sse2
    friend mdp_matrix operator/(const mdp_matrix &a, const mdp_matrix &b);

    friend mdp_matrix operator+(const mdp_matrix &a, mdp_complex b);
    friend mdp_matrix operator-(const mdp_matrix &a, mdp_complex b);
    friend mdp_matrix operator*(const mdp_matrix &a, mdp_complex b);
    friend mdp_matrix operator/(const mdp_matrix &a, mdp_complex b);

    friend mdp_matrix operator+(mdp_complex a, const mdp_matrix &b);
    friend mdp_matrix operator-(mdp_complex a, const mdp_matrix &b);
    friend mdp_matrix operator*(mdp_complex a, const mdp_matrix &b);
    friend mdp_matrix operator/(mdp_complex a, const mdp_matrix &b);

    friend mdp_matrix operator+(const mdp_matrix &a, mdp_real b);
    friend mdp_matrix operator-(const mdp_matrix &a, mdp_real b);
    friend mdp_matrix operator*(const mdp_matrix &a, mdp_real b);
    friend mdp_matrix operator/(const mdp_matrix &a, mdp_real b);

    friend mdp_matrix operator+(mdp_real a, const mdp_matrix &b);
    friend mdp_matrix operator-(mdp_real a, const mdp_matrix &b);
    friend mdp_matrix operator*(mdp_real a, const mdp_matrix &b);
    friend mdp_matrix operator/(mdp_real a, const mdp_matrix &b);

    friend mdp_matrix inv(const mdp_matrix &a);
    friend mdp_matrix pow(const mdp_matrix &a, mdp_uint b);
    friend mdp_matrix exp(const mdp_matrix &a);
    friend mdp_matrix log(const mdp_matrix &a);
    friend mdp_matrix sin(const mdp_matrix &a);
    friend mdp_matrix cos(const mdp_matrix &a);

    friend mdp_matrix mdp_identity();
    friend mdp_matrix mdp_zero();

    friend mdp_real max(const mdp_matrix &a);
    friend mdp_matrix submatrix(const mdp_matrix &a, mdp_uint i, mdp_uint j);
    friend mdp_complex det(const mdp_matrix &a);
    friend mdp_complex trace(const mdp_matrix &a);
    friend mdp_matrix hermitian(const mdp_matrix &a);
    friend mdp_matrix transpose(const mdp_matrix &a);
    friend mdp_matrix conj(const mdp_matrix &a);
  };

  std::ostream &operator<<(std::ostream &os, const mdp_matrix &a)
  {
    for (mdp_uint i = 0; i < a.rows(); i++)
    {
      if (i == 0)
        os << "[[";
      else
        os << " [";
      for (mdp_uint j = 0; j < a.cols(); j++)
        if (j == 0)
          os << " " << a(i, j);
        else
          os << ", " << a(i, j) << " ";
      if (i == (a.rows() - 1))
        os << "]]\n";
      else
        os << "], \n";
    }
    return os;
  }

  mdp_matrix operator+(const mdp_matrix &a)
  {
    return a;
  }

  mdp_matrix operator-(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.rows(), a.cols());

    for (mdp_uint i = 0; i < a.rows(); i++)
      for (mdp_uint j = 0; j < a.cols(); j++)
        tmp(i, j) = -a(i, j);
    return tmp;
  }

  mdp_matrix operator+(const mdp_matrix &x, const mdp_matrix &y)
  {
    mdp_matrix z(x.rows(), x.cols());
#ifdef CHECK_ALL
    if (x.rows() != y.rows() || x.cols() != y.cols())
      error("mdp_matrix::operator+()\nWrong argument size");
#endif
    for (mdp_uint i = 0; i < z.size(); i++)
      z[i] = x[i] + y[i];
    return z;
  }

  mdp_matrix operator-(const mdp_matrix &x, const mdp_matrix &y)
  {
    if (x.rows() != y.rows() || x.cols() != y.cols())
      error("mdp_matrix::operator-()\nwrong argument size");
    mdp_matrix z(x.rows(), x.cols());
    for (mdp_uint i = 0; i < z.size(); i++)
      z[i] = x[i] - y[i];
    return z;
  }

  mdp_matrix operator*(const mdp_matrix &x, const mdp_matrix &y)
  {
    if (x.cols() != y.rows())
      error("mdp_matrix::operator*()\nwrong argument size");
    mdp_matrix z(x.rows(), y.cols());
#if defined(MATRIX_SSE2) && defined(USE_DOUBLE_PRECISION)
    if (x.rows() == x.cols() && y.rows() == 3)
    {
      _sse_su3 *u = (_sse_su3 *)x.m_data;
      _sse_double *in = (_sse_double *)y.m_data;
      _sse_double *out = (_sse_double *)z.m_data;
      for (mdp_uint i = 0; i < x.cols(); i++, in++, out++)
      {
        _sse_double_load_123(*in, *(in + x.cols()), *(in + 2 * x.cols()));
        _sse_double_su3_multiply(*u);
        _sse_double_store_up_123(*out, *(out + x.cols()), *(out + 2 * x.cols()));
      }
      return z;
    }
#endif
#if defined(MATRIX_SSE2) && !defined(USE_DOUBLE_PRECISION)
    if (x.rows() == x.cols() && y.size() == 3)
    {
      _sse_su3 *u = (_sse_su3 *)x.m_data;
      _sse_su3_vector *in = (_sse_su3_vector *)y.m_data;
      _sse_su3_vector *out = (_sse_su3_vector *)z.m_data;
      _sse_float_pair_load(*in, *in);
      _sse_float_su3_multiply(*u);
      _sse_float_pair_store_up(*out, *out);
      return z;
    }
#endif
    for (mdp_uint i = 0; i < x.rows(); i++)
      for (mdp_uint j = 0; j < y.cols(); j++)
      {
        z(i, j) = x(i, 0) * y(0, j);
        for (mdp_uint k = 1; k < x.cols(); k++)
          z(i, j) += x(i, k) * y(k, j);
      }
    return z;
  }

  mdp_matrix operator/(const mdp_matrix &a, const mdp_matrix &b)
  {
    return a * inv(b);
  }

  mdp_matrix operator+(const mdp_matrix &a, mdp_complex b)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator+(...)\nmdp_matrix is not squared");
    mdp_matrix tmp;
    tmp = a;
    for (mdp_uint i = 0; i < a.cols(); i++)
      tmp(i, i) += b;
    return tmp;
  }

  mdp_matrix operator-(const mdp_matrix &a, mdp_complex b)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
    mdp_matrix tmp;

    tmp = a;
    for (mdp_uint i = 0; i < a.cols(); i++)
      tmp(i, i) -= b;
    return tmp;
  }

  mdp_matrix operator*(const mdp_matrix &y, mdp_complex x)
  {
    mdp_matrix z(y.rows(), y.cols());
#ifdef MATRIX_SSE2
    if (y.rows() == 3)
    {
      static _sse_float factor1 ALIGN16;
      static _sse_float factor2 ALIGN16;
      static _sse_double factor3 ALIGN16;
      static _sse_double factor4 ALIGN16;
      _sse_su3_vector *in = (_sse_su3_vector *)y.m_data;
      _sse_su3_vector *out = (_sse_su3_vector *)z.m_data;
#ifdef USE_DOUBLE_PRECISION
      factor3.c1 = factor3.c2 = x.imag();
      factor4.c1 = factor4.c2 = x.real() / x.imag();
      for (mdp_uint i = 0; i < y.cols(); i++, in++, out++)
      {
        _sse_double_load(*in);
        _sse_double_vector_mulc(factor3, factor4);
        _sse_double_store(*out);
      }
#else
      factor1.c1 = factor1.c2 = factor1.c3 = factor1.c4 = x.imag();
      factor2.c1 = factor2.c2 = factor2.c3 = factor2.c4 = x.real() / x.imag();
      mdp_uint i = 0;
      for (i = 0; i < y.cols() - 1; i += 2, in += 2, out += 2)
      {
        _sse_float_pair_load(*in, *(in + 1));
        _sse_float_vector_mulc(factor1, factor2);
        _sse_float_pair_store(*out, *(out + 1));
      }
      if (i == y.cols() - 1)
      {
        _sse_float_pair_load(*in, *in);
        _sse_float_vector_mulc(factor1, factor2);
        _sse_float_pair_store(*out, *out);
      }
#endif
      return z;
    }
#endif
    for (mdp_uint i = 0; i < y.size(); i++)
      z[i] = x * y[i];
    return z;
  }

  mdp_matrix operator/(const mdp_matrix &a, mdp_complex b)
  {
    return a * (1.0 / b);
  }

  mdp_matrix operator+(mdp_complex b, const mdp_matrix &a)
  {
    return a + b;
  }

  mdp_matrix operator-(mdp_complex b, const mdp_matrix &a)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
    mdp_matrix tmp(a.rows(), a.cols());

    for (mdp_uint i = 0; i < a.rows(); i++)
    {
      for (mdp_uint j = 0; j < a.cols(); j++)
        tmp(i, j) = -a(i, j);
      tmp(i, i) += b;
    }
    return tmp;
  }

  mdp_matrix operator*(mdp_complex x, const mdp_matrix &y)
  {
    return y * x;
  }

  mdp_matrix operator/(mdp_complex b, const mdp_matrix &a)
  {
    return b * inv(a);
  }

  mdp_matrix operator+(const mdp_matrix &a, mdp_real b)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator+(...)\nmdp_matrix is not squared");
    mdp_matrix tmp;

    tmp = a;
    for (mdp_uint i = 0; i < a.cols(); i++)
      tmp(i, i).real() += b;
    return tmp;
  }

  mdp_matrix operator-(const mdp_matrix &a, mdp_real b)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
    mdp_matrix tmp;

    tmp = a;
    for (mdp_uint i = 0; i < a.cols(); i++)
      tmp(i, i).real() -= b;
    return tmp;
  }

  mdp_matrix operator*(const mdp_matrix &y, mdp_real x)
  {
    mdp_matrix z(y.rows(), y.cols());
#ifdef MATRIX_SSE2
    if (y.rows() == 3)
    {
      static _sse_float factor1 ALIGN16;
      static _sse_double factor2 ALIGN16;

      _sse_su3_vector *in = (_sse_su3_vector *)y.m_data;
      _sse_su3_vector *out = (_sse_su3_vector *)z.m_data;
#ifdef USE_DOUBLE_PRECISION
      factor2.c1 = factor2.c2 = x;
      for (mdp_uint i = 0; i < y.cols(); i++, in++, out++)
      {
        _sse_double_load(*in);
        _sse_double_vector_mul(factor2);
        _sse_double_store(*out);
      }
#else
      factor1.c1 = factor1.c2 = factor1.c3 = factor1.c4 = x;
      mdp_uint i = 0;
      for (i = 0; i < y.cols() - 1; i += 2, in += 2, out += 2)
      {
        _sse_float_pair_load(*in, *(in + 1));
        _sse_float_vector_mul(factor1);
        _sse_float_pair_store(*out, *(out + 1));
      }
      if (i == y.cols() - 1)
      {
        _sse_float_pair_load(*in, *in);
        _sse_float_vector_mul(factor1);
        _sse_float_pair_store(*out, *out);
      }
#endif
      return z;
    }
#endif
    for (mdp_uint i = 0; i < y.size(); i++)
      z[i] = x * y[i];
    return z;
  }

  mdp_matrix operator/(const mdp_matrix &a, mdp_real b)
  {
    return a * (1.0 / b);
  }

  mdp_matrix operator+(mdp_real b, const mdp_matrix &a)
  {
    return a + b;
  }

  mdp_matrix operator-(mdp_real b, const mdp_matrix &a)
  {
    if (a.cols() != a.rows())
      error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
    mdp_matrix tmp(a.rows(), a.cols());

    for (mdp_uint i = 0; i < a.rows(); i++)
    {
      for (mdp_uint j = 0; j < a.cols(); j++)
        tmp(i, j) = -a(i, j);
      tmp(i, i) += b;
    }
    return tmp;
  }

  mdp_matrix operator*(mdp_real a, const mdp_matrix &b)
  {
    return b * a;
  }

  mdp_matrix operator/(mdp_real b, const mdp_matrix &a)
  {
    return b * inv(a);
  }

  /** @brief Create square identity matrix of size \e i
   *
   * @param i size of the matrix
   *
   * @return Identity matrix
   */
  mdp_matrix mdp_identity(mdp_uint i)
  {
    mdp_matrix tmp(i, i);
    for (mdp_uint r = 0; r < i; r++)
      for (mdp_uint c = 0; c < i; c++)
        tmp(r, c) = (r == c) ? mdp_complex(1, 0) : mdp_complex(0, 0);
    return tmp;
  }

  /** @brief Create a matrix of size \e ixj filled with zeros
   *
   * @param i number of rows
   * @param j number of columns
   *
   * @return Matrix filled with zeros
   */
  mdp_matrix mdp_zero(mdp_uint i, mdp_uint j)
  {
    mdp_matrix tmp(i, j);
    for (mdp_uint r = 0; r < i; r++)
      for (mdp_uint c = 0; c < j; c++)
        tmp(r, c) = mdp_complex(0, 0);
    return tmp;
  }

  /** @brief Create square matrix of size \e i filled with zeros
   *
   * @param i size of the matrix
   *
   * @return Matrix filled with zeros
   */
  mdp_matrix mdp_zero(mdp_uint i)
  {
    return mdp_zero(i, i);
  }

  /** @brief Return matrix element with the highest value
   *
   * @param a Complex matrix
   *
   * @return Element a(i,j) with the highest absolute value
   */
  mdp_real max(const mdp_matrix &a)
  {
    double x = 0, y = 0;
    for (mdp_uint i = 0; i < a.size(); i++)
      if ((y = abs(a[i])) > x)
        x = y;
    return x;
  }

  /** @brief Return matrix without i-th row and j-th column
   *
   * @param a Input matrix
   *
   * @param i Row to be excluded from the matrix
   *
   * @param j Column to be excluded from the matrix
   *
   * @return Matrix without i-th row and j-th column
   */
  mdp_matrix submatrix(const mdp_matrix &a, mdp_uint i, mdp_uint j)
  {
#ifdef CHECK_ALL
    if (((a.rows() < 2) || (a.cols() < 2)) ||
        ((a.rows() - 1 < i) || (a.cols() - 1 < j)))
      error("submatrix(...)\nWrong dimensions in submatrix");
#endif
    mdp_matrix tmp(a.rows() - 1, a.cols() - 1);
    for (mdp_uint r = 0; r < i; r++)
      for (mdp_uint c = 0; c < j; c++)
        tmp(r, c) = a(r, c);
    for (mdp_uint r = i + 1; r < a.rows(); r++)
      for (mdp_uint c = 0; c < j; c++)
        tmp(r - 1, c) = a(r, c);
    for (mdp_uint r = 0; r < i; r++)
      for (mdp_uint c = j + 1; c < a.cols(); c++)
        tmp(r, c - 1) = a(r, c);
    for (mdp_uint r = i + 1; r < a.rows(); r++)
      for (mdp_uint c = j + 1; c < a.cols(); c++)
        tmp(r - 1, c - 1) = a(r, c);
    return tmp;
  }

  mdp_complex det(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("det(...)\nmdp_matrix is not squared");
#endif
    if (a.rows() == 0)
      return 0;
    if (a.rows() == 1)
      return a(0, 0);
    mdp_uint i, j, k, l;
    mdp_matrix A;
    A = a;
    mdp_uint cols = a.cols();
    mdp_uint rows = a.rows();
    mdp_complex tmp, pivot, x = mdp_complex(1, 0);
    for (i = 0; i < cols; i++)
    {
      for (j = i; (A(i, j) == mdp_complex(0, 0)) && (j < cols); j++)
        ;
      if (j == cols)
        return 0;
      if (i != j)
      {
        for (k = 0; k < rows; k++)
        {
          tmp = A(k, j);
          A(k, j) = A(k, i);
          A(k, i) = tmp;
        }
        x *= -A(i, i);
      }
      else
        x *= A(i, i);
      for (k = i + 1; k < rows; k++)
      {
        pivot = A(k, i) / A(i, i);
        for (l = i; l < cols; l++)
          A(k, l) -= pivot * A(i, l);
      }
    }
    return x;
  }

  mdp_matrix inv(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if ((a.rows() != a.cols()) || (a.rows() == 0))
      error("inv(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tma, tmp;
    mdp_complex x, pivot;
    mdp_uint rmax;
    tma = a;
    tmp = mdp_identity(a.rows());
    for (mdp_uint c = 0; c < a.cols(); c++)
    {
      rmax = c;
      pivot = tma(c, c);
      for (mdp_uint r = c + 1; r < a.rows(); r++)
        if (abs(tma(r, c)) > abs(pivot))
        {
          rmax = r;
          pivot = tma(r, c);
        }
      for (mdp_uint i = 0; i < a.cols(); i++)
      {
        x = tma(rmax, i);
        tma(rmax, i) = tma(c, i);
        tma(c, i) = x / pivot;
        x = tmp(rmax, i);
        tmp(rmax, i) = tmp(c, i);
        tmp(c, i) = x / pivot;
      }
      for (mdp_uint r = 0; r < a.rows(); r++)
        if (r != c)
        {
          pivot = tma(r, c);
          for (mdp_uint i = 0; i < a.cols(); i++)
          {
            tma(r, i) -= pivot * tma(c, i);
            tmp(r, i) -= pivot * tmp(c, i);
          }
        }
    }
    return tmp;
  }

  mdp_matrix pow(const mdp_matrix &a, int i)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("pow(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tmp;
    tmp = mdp_identity(a.cols());
    mdp_int j = (i < 0) ? -i : i;
    for (; j > 0; j--)
      tmp = tmp * a;
    if (i < 0)
      tmp = inv(tmp);
    return tmp;
  }

  mdp_matrix exp(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("exp(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tmp;
    mdp_matrix term;
    mdp_uint i = 1;
    term = a;
    tmp = mdp_identity(a.rows());
    tmp += a;
    do
    {
      term = (1. / ++i) * term * a;
      tmp += term;
    } while (max(term) > mdp_precision);
    return tmp;
  }

  mdp_matrix log(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("log(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tmp, b, c, t1;
    mdp_int i = 1;
    mdp_uint j = 1;
    b = mdp_identity(a.cols());
    b = a - b;
    c = b;
    tmp = b;
    do
    {
      c = c * b;
      t1 = ((mdp_real)(i = -i) / (j += 1)) * c;
      tmp += t1;
    } while (max(t1) > mdp_precision);
    return tmp;
  }

  mdp_matrix sin(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("sin(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tmp, t1;
    mdp_uint i = 1;
    t1 = a;
    tmp = t1;
    // pruintf("\n");
    do
    {
      t1 = ((mdp_real)-1.0 / (++i)) * t1 * a * a;
      t1 *= (mdp_real)1.0 / (++i);
      tmp += t1;
    } while (max(t1) > mdp_precision);
    return tmp;
  }

  mdp_matrix cos(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("cos(...)\nmdp_matrix is not squared");
#endif
    mdp_matrix tmp, t1;
    mdp_uint i = 0;
    t1 = mdp_identity(a.rows());
    tmp = t1;
    do
    {
      t1 = ((mdp_real)-1.0 / (++i)) * t1 * a * a;
      t1 *= (mdp_real)1.0 / (++i);
      tmp += t1;
    } while (max(t1) > mdp_precision);
    return tmp;
  }

  mdp_complex trace(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("trace(...)\nmdp_matrix is not squared");
#endif
    mdp_complex x;
    for (mdp_uint c = 0; c < a.cols(); c++)
      x += a(c, c);
    return x;
  }

  mdp_matrix transpose(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.cols(), a.rows());

    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
        tmp(c, r) = a(r, c);
    return tmp;
  }

  mdp_matrix hermitian(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.cols(), a.rows());

#if defined(MATRIX_SSE2) && defined(USE_DOUBLE_PRECISION)
    if (a.cols() == 3 && a.rows() == 3)
    {
      _sse_double_hermitian_su3((_sse_double *)tmp.address(),
                                (_sse_double *)a.address());
      return tmp;
    }
#endif

    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
        tmp(c, r) = conj(a(r, c));
    return tmp;
  }

  mdp_matrix conj(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.rows(), a.cols());
    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
        tmp(r, c) = conj(a(r, c));
    return tmp;
  }
} // namespace MDP

#endif /* MDP_MATRIX_ */
