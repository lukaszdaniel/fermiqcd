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
#include <span>
#include <algorithm>
#include <type_traits>
#include "mdp_macros.h"
#include "mdp_global_vars.h"
#include "mdp_communicator.h"
#include "mdp_complex.h"

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
    mdp_uint m_rows;
    mdp_uint m_cols;
    std::unique_ptr<mdp_complex[]> m_owner;
    std::span<mdp_complex> m_data;

    void reallocate()
    {
      if (size() == 0)
      {
        m_owner.reset();
        m_data = {};
        return;
      }

      m_owner = std::make_unique<mdp_complex[]>(size());
      m_data = std::span<mdp_complex>(m_owner.get(), size());

      std::fill(m_data.begin(), m_data.end(), mdp_complex{});
    }

  public:
    /** @brief Create a matrix
     *
     * @param r Number of rows. Default is 3.
     * @param c Number of columns. Default is 3.
     */
    mdp_matrix(mdp_uint r = 3, mdp_uint c = 3)
        : m_rows(r), m_cols(c),
          m_owner(std::make_unique<mdp_complex[]>(r * c)),
          m_data(m_owner.get(), r * c)
    {
      std::fill(m_data.begin(), m_data.end(), mdp_complex{});
    }

    /** @brief Copy constructor
     *
     * @param a Matrix to be copied.
     */
    mdp_matrix(const mdp_matrix &a)
        : m_rows(a.m_rows),
          m_cols(a.m_cols),
          m_owner(std::make_unique<mdp_complex[]>(a.size())),
          m_data(m_owner.get(), a.size())
    {
      std::copy(a.m_data.begin(), a.m_data.end(), m_data.begin());
    }

    /** @brief View constructor used to format array of complex numbers
     *
     * @param z Pointer to the array of complex number
     * @param r Number of rows array needs to be split into.
     * @param c Number of columns array needs to be split into.
     *
     * @note Length of \e z array needs to be equal to \e rc.
     */
    mdp_matrix(mdp_complex *ptr, mdp_uint r, mdp_uint c)
        : m_rows(r), m_cols(c),
          m_owner(nullptr),
          m_data(ptr, r * c)
    {
    }

    ~mdp_matrix() = default;

    mdp_matrix(mdp_matrix &&) noexcept = default;
    // mdp_matrix &operator=(mdp_matrix &&) noexcept = default;

    /** @brief Set matrix size
     *
     * @param r Number of rows for the matrix.
     * @param c Number of columns for the matrix.
     */
    void dimension(mdp_uint r, mdp_uint c)
    {
      m_rows = r;
      m_cols = c;
      reallocate();
    }

    mdp_matrix &operator=(const mdp_matrix &x)
    {
      if (this != &x)
      {
        if (m_rows != x.m_rows || m_cols != x.m_cols)
        {
          m_rows = x.m_rows;
          m_cols = x.m_cols;
          reallocate();
        }
        std::copy(x.m_data.begin(), x.m_data.end(), m_data.begin());
      }

      return *this;
    }

    auto begin() { return m_data.begin(); }
    auto end() { return m_data.end(); }
    auto begin() const { return m_data.begin(); }
    auto end() const { return m_data.end(); }
    mdp_complex *data() { return m_data.data(); }
    const mdp_complex *data() const { return m_data.data(); }

    /** @brief returns the (i,j) element of the matrix
     */
    inline mdp_complex &operator()(mdp_uint i, mdp_uint j)
    {
      return m_data[i * m_cols + j];
    }

    /** @brief returns the (i,j) const element of the matrix
     */
    inline const mdp_complex &operator()(mdp_uint i, mdp_uint j) const
    {
      return m_data[i * m_cols + j];
    }

    /** @brief returns the i-th row of the matrix aligned as column vector
     */
    mdp_matrix operator()(mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= m_rows)
        error("mdp_matrix::operator()\nIndex out of bounds");
#endif
      return mdp_matrix(m_data.data() + i * m_cols, m_cols, 1);
    }

    /** @brief Begining of this matrix
     *
     * @return Pointer to the begining of this matrix
     */
    mdp_complex *address() const
    {
      return m_data.data();
    }

    mdp_complex &operator[](mdp_uint i)
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_matrix::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    const mdp_complex &operator[](mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_matrix::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    /** @brief The number of rows
     *
     * @return The number of rows
     */
    mdp_uint rows() const
    {
      return m_rows;
    }

    /** @brief The number of columns
     *
     * @return The number of columns
     */
    mdp_uint cols() const
    {
      return m_cols;
    }

    /** @brief Total size of this matrix
     *
     * @return Total number of elements
     */
    mdp_uint size() const
    {
      return m_rows * m_cols;
    }

    mdp_matrix operator+() const
    {
      return (*this);
    }

    mdp_matrix operator-() const
    {
      mdp_matrix tmp(m_rows, m_cols);

      for (mdp_uint i = 0; i < size(); i++)
        tmp[i] = -m_data[i];

      return tmp;
    }

    mdp_matrix operator+(const mdp_matrix &x) const
    {
#ifdef CHECK_ALL
      if (m_rows != x.m_rows || m_cols != x.m_cols)
        error("mdp_matrix::operator+()\nWrong argument size");
#endif
      mdp_matrix z(m_rows, m_cols);
      for (mdp_uint i = 0; i < size(); i++)
      {
        z[i] = m_data[i] + x.m_data[i];
      }
      return z;
    }

    mdp_matrix operator-(const mdp_matrix &x) const
    {
#ifdef CHECK_ALL
      if (m_rows != x.m_rows || m_cols != x.m_cols)
        error("mdp_matrix::operator-()\nWrong argument size");
#endif
      mdp_matrix z(m_rows, m_cols);
      for (mdp_uint i = 0; i < size(); i++)
      {
        z[i] = m_data[i] - x.m_data[i];
      }
      return z;
    }

    mdp_matrix operator*(const mdp_matrix &x) const
    {
#ifdef CHECK_ALL
      if (m_cols != x.m_rows)
        error("mdp_matrix::operator*()\nWrong argument size");
#endif

      if (m_rows == 3 && m_cols == 3 && x.m_rows == 3 && x.m_cols == 3) [[likely]]
      {
        // optimized multiplication for 3x3 matrices
        mdp_matrix C(3, 3);

        const mdp_complex *A = m_data.data();
        const mdp_complex *b = x.m_data.data();
        mdp_complex *c = C.m_data.data();

        const mdp_complex &a00 = A[0];
        const mdp_complex &a01 = A[1];
        const mdp_complex &a02 = A[2];
        const mdp_complex &a10 = A[3];
        const mdp_complex &a11 = A[4];
        const mdp_complex &a12 = A[5];
        const mdp_complex &a20 = A[6];
        const mdp_complex &a21 = A[7];
        const mdp_complex &a22 = A[8];

        const mdp_complex &b00 = b[0];
        const mdp_complex &b01 = b[1];
        const mdp_complex &b02 = b[2];
        const mdp_complex &b10 = b[3];
        const mdp_complex &b11 = b[4];
        const mdp_complex &b12 = b[5];
        const mdp_complex &b20 = b[6];
        const mdp_complex &b21 = b[7];
        const mdp_complex &b22 = b[8];

        c[0] = a00 * b00 + a01 * b10 + a02 * b20;
        c[1] = a00 * b01 + a01 * b11 + a02 * b21;
        c[2] = a00 * b02 + a01 * b12 + a02 * b22;

        c[3] = a10 * b00 + a11 * b10 + a12 * b20;
        c[4] = a10 * b01 + a11 * b11 + a12 * b21;
        c[5] = a10 * b02 + a11 * b12 + a12 * b22;

        c[6] = a20 * b00 + a21 * b10 + a22 * b20;
        c[7] = a20 * b01 + a21 * b11 + a22 * b21;
        c[8] = a20 * b02 + a21 * b12 + a22 * b22;

        return C;
      }

      mdp_matrix z(m_rows, x.m_cols);
      for (mdp_uint i = 0; i < m_rows; i++)
      {
        for (mdp_uint k = 0; k < m_cols; k++)
        {
          const mdp_complex aik = m_data[i * m_cols + k];

          for (mdp_uint j = 0; j < m_cols; j++)
          {
            z.m_data[i * m_cols + j] += aik * x.m_data[k * m_cols + j];
          }
        }
      }
      return z;
    }

    mdp_matrix operator/(const mdp_matrix &a) const
    {
      return (*this) * a.inv();
    }

    mdp_matrix &operator+=(const mdp_matrix &a)
    {
#ifdef CHECK_ALL
      if (m_rows != a.m_rows || m_cols != a.m_cols)
        error("mdp_matrix::operator+=()\nWrong argument size");
#endif

      for (mdp_uint i = 0; i < size(); ++i)
        m_data[i] += a.m_data[i];

      return *this;
    }

    mdp_matrix &operator-=(const mdp_matrix &a)
    {
#ifdef CHECK_ALL
      if (m_rows != a.m_rows || m_cols != a.m_cols)
        error("mdp_matrix::operator-=()\nWrong argument size");
#endif

      for (mdp_uint i = 0; i < size(); ++i)
        m_data[i] -= a.m_data[i];
      return *this;
    }

    mdp_matrix &operator*=(const mdp_matrix &a)
    {
      *this = (*this) * a;
      return *this;
    }

    mdp_matrix &operator/=(const mdp_matrix &a)
    {
      *this = (*this) / a;
      return *this;
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix operator+(T b) const
    {
      if (m_cols != m_rows)
        error("mdp_matrix::operator+(...)\nmdp_matrix is not squared");

      mdp_matrix tmp((*this));

      for (mdp_uint i = 0; i < m_rows; i++)
        tmp(i, i) += b;

      return tmp;
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix operator-(T b) const
    {
      if (m_cols != m_rows)
        error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");

      mdp_matrix tmp((*this));

      for (mdp_uint i = 0; i < m_rows; i++)
        tmp(i, i) -= b;

      return tmp;
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix operator*(T x) const
    {
      mdp_matrix z(m_rows, m_cols);

      // // optimized multiplication for 3x3 matrix
      if (m_rows == 3 && m_cols == 3) [[likely]]
      {
        const mdp_complex &a00 = m_data[0];
        const mdp_complex &a01 = m_data[1];
        const mdp_complex &a02 = m_data[2];
        const mdp_complex &a10 = m_data[3];
        const mdp_complex &a11 = m_data[4];
        const mdp_complex &a12 = m_data[5];
        const mdp_complex &a20 = m_data[6];
        const mdp_complex &a21 = m_data[7];
        const mdp_complex &a22 = m_data[8];

        z.m_data[0] = a00 * x;
        z.m_data[1] = a01 * x;
        z.m_data[2] = a02 * x;
        z.m_data[3] = a10 * x;
        z.m_data[4] = a11 * x;
        z.m_data[5] = a12 * x;
        z.m_data[6] = a20 * x;
        z.m_data[7] = a21 * x;
        z.m_data[8] = a22 * x;

        return z;
      }

      for (mdp_uint i = 0; i < size(); i++)
        z[i] = x * m_data[i];

      return z;
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix operator/(T b) const
    {
      return (*this) * (mdp_complex(1) / mdp_complex(b));
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix &operator+=(T x)
    {
      mdp_complex c = mdp_complex(x);

      for (mdp_uint i = 0; i < m_rows; ++i)
        (*this)(i, i) += c;

      return *this;
    }

    template <class T>
      requires std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>
    mdp_matrix &operator-=(T x)
    {
      mdp_complex c = mdp_complex(x);

      for (mdp_uint i = 0; i < m_rows; ++i)
        (*this)(i, i) -= c;

      return *this;
    }

    template <class T>
      requires(std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>)
    mdp_matrix &operator*=(T x)
    {
      mdp_complex c = mdp_complex(x);

      for (mdp_uint i = 0; i < size(); ++i)
        m_data[i] *= c;

      return *this;
    }

    template <class T>
      requires(std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>)
    mdp_matrix &operator/=(T x)
    {
      mdp_complex c = mdp_complex(x);

      (*this) = (*this) / c;
      return *this;
    }

    template <class T>
      requires(std::is_arithmetic_v<T> || std::same_as<T, mdp_complex>)
    mdp_matrix &operator=(T x)
    {
      mdp_complex c = mdp_complex(x);

      for (mdp_uint i = 0; i < rows(); i++)
        for (mdp_uint j = 0; j < cols(); j++)
        {
          (*this)(i, j) = (i == j) ? c : mdp_complex(0);
        }
      return *this;
    }

    mdp_complex det() const
    {
#ifdef CHECK_ALL
      if (m_rows != m_cols)
        error("det(...)\nmdp_matrix is not squared");
#endif

      const mdp_uint n = m_rows;

      // Handle most common cases quickly here.
      if (n == 0)
        return mdp_complex(0);

      if (n == 1)
        return m_data[0];

      if (n == 2)
      {
        const mdp_complex a = m_data[0];
        const mdp_complex b = m_data[1];
        const mdp_complex c = m_data[2];
        const mdp_complex d = m_data[3];

        return a * d - b * c;
      }

      if (n == 3) [[likely]]
      {
        const mdp_complex a00 = m_data[0];
        const mdp_complex a01 = m_data[1];
        const mdp_complex a02 = m_data[2];

        const mdp_complex a10 = m_data[3];
        const mdp_complex a11 = m_data[4];
        const mdp_complex a12 = m_data[5];

        const mdp_complex a20 = m_data[6];
        const mdp_complex a21 = m_data[7];
        const mdp_complex a22 = m_data[8];

        return a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20);
      }

      if (n == 4)
      {
        const auto &a = m_data;

        mdp_complex a0 = a[0] * a[5] - a[1] * a[4];
        mdp_complex a1 = a[0] * a[6] - a[2] * a[4];
        mdp_complex a2 = a[0] * a[7] - a[3] * a[4];
        mdp_complex a3 = a[1] * a[6] - a[2] * a[5];
        mdp_complex a4 = a[1] * a[7] - a[3] * a[5];
        mdp_complex a5 = a[2] * a[7] - a[3] * a[6];

        mdp_complex b0 = a[8] * a[13] - a[9] * a[12];
        mdp_complex b1 = a[8] * a[14] - a[10] * a[12];
        mdp_complex b2 = a[8] * a[15] - a[11] * a[12];
        mdp_complex b3 = a[9] * a[14] - a[10] * a[13];
        mdp_complex b4 = a[9] * a[15] - a[11] * a[13];
        mdp_complex b5 = a[10] * a[15] - a[11] * a[14];

        return a0 * b5 - a1 * b4 + a2 * b3 +
               a3 * b2 - a4 * b1 + a5 * b0;
      }

      mdp_matrix A(*this);

      mdp_complex det = mdp_complex(1);
      int sign = 1;

      for (mdp_uint k = 0; k < n; ++k)
      {
        mdp_uint pivot_row = k;
        auto pivot_norm = abs2(A(k, k));

        for (mdp_uint r = k + 1; r < n; ++r)
        {
          auto v = abs2(A(r, k));
          if (v > pivot_norm)
          {
            pivot_norm = v;
            pivot_row = r;
          }
        }

        if (pivot_norm == decltype(pivot_norm)(0))
          return mdp_complex(0);

        if (pivot_row != k)
        {
          std::swap_ranges(&A(pivot_row, 0), &A(pivot_row, 0) + n, &A(k, 0));
          sign = -sign;
        }

        mdp_complex pivot = A(k, k);
        det *= pivot;

        mdp_complex invPivot = mdp_complex(1) / pivot;

        for (mdp_uint i = k + 1; i < n; ++i)
        {
          mdp_complex factor = A(i, k) * invPivot;

          for (mdp_uint j = k + 1; j < n; ++j)
            A(i, j) -= factor * A(k, j);
        }
      }

      return sign * det;
    }

    mdp_matrix inv() const
    {
#ifdef CHECK_ALL
      if (m_rows != m_cols)
        error("inv(...)\nmdp_matrix is not squared");
      if (m_rows == 0)
        error("inv(...)\nmdp_matrix is of size 0");
#endif

      const mdp_uint n = m_rows;
      mdp_matrix ans(n, n);

      // Handle most common cases quickly here.
      if (n == 1)
      {
        if (m_data[0] == mdp_complex(0)) [[unlikely]]
          error("inv(...)\ndeterminant is zero");

        ans.m_data[0] = mdp_complex(1) / m_data[0];
        return ans;
      }

      if (n == 2)
      {
        const mdp_complex a = m_data[0];
        const mdp_complex b = m_data[1];
        const mdp_complex c = m_data[2];
        const mdp_complex d = m_data[3];

        const mdp_complex det = a * d - b * c;

        if (det == mdp_complex(0)) [[unlikely]]
          error("inv(...)\ndeterminant is zero");

        const mdp_complex invdet = mdp_complex(1) / det;

        ans.m_data[0] = d * invdet;
        ans.m_data[1] = -b * invdet;
        ans.m_data[2] = -c * invdet;
        ans.m_data[3] = a * invdet;

        return ans;
      }

      if (n == 3) [[likely]]
      {
        const mdp_complex a00 = m_data[0];
        const mdp_complex a01 = m_data[1];
        const mdp_complex a02 = m_data[2];

        const mdp_complex a10 = m_data[3];
        const mdp_complex a11 = m_data[4];
        const mdp_complex a12 = m_data[5];

        const mdp_complex a20 = m_data[6];
        const mdp_complex a21 = m_data[7];
        const mdp_complex a22 = m_data[8];

        const mdp_complex c00 = a11 * a22 - a12 * a21;
        const mdp_complex c01 = a02 * a21 - a01 * a22;
        const mdp_complex c02 = a01 * a12 - a02 * a11;

        const mdp_complex c10 = a12 * a20 - a10 * a22;
        const mdp_complex c11 = a00 * a22 - a02 * a20;
        const mdp_complex c12 = a02 * a10 - a00 * a12;

        const mdp_complex c20 = a10 * a21 - a11 * a20;
        const mdp_complex c21 = a01 * a20 - a00 * a21;
        const mdp_complex c22 = a00 * a11 - a01 * a10;

        const mdp_complex det = a00 * c00 + a01 * c10 + a02 * c20;

        if (det == mdp_complex(0)) [[unlikely]]
          error("inv(...)\ndeterminant is zero");

        const mdp_complex invdet = mdp_complex(1) / det;

        ans.m_data[0] = c00 * invdet;
        ans.m_data[1] = c01 * invdet;
        ans.m_data[2] = c02 * invdet;

        ans.m_data[3] = c10 * invdet;
        ans.m_data[4] = c11 * invdet;
        ans.m_data[5] = c12 * invdet;

        ans.m_data[6] = c20 * invdet;
        ans.m_data[7] = c21 * invdet;
        ans.m_data[8] = c22 * invdet;

        return ans;
      }

      if (n == 4)
      {
        const auto &a = m_data;

        mdp_complex a0 = a[0] * a[5] - a[1] * a[4];
        mdp_complex a1 = a[0] * a[6] - a[2] * a[4];
        mdp_complex a2 = a[0] * a[7] - a[3] * a[4];
        mdp_complex a3 = a[1] * a[6] - a[2] * a[5];
        mdp_complex a4 = a[1] * a[7] - a[3] * a[5];
        mdp_complex a5 = a[2] * a[7] - a[3] * a[6];

        mdp_complex b0 = a[8] * a[13] - a[9] * a[12];
        mdp_complex b1 = a[8] * a[14] - a[10] * a[12];
        mdp_complex b2 = a[8] * a[15] - a[11] * a[12];
        mdp_complex b3 = a[9] * a[14] - a[10] * a[13];
        mdp_complex b4 = a[9] * a[15] - a[11] * a[13];
        mdp_complex b5 = a[10] * a[15] - a[11] * a[14];

        mdp_complex det =
            a0 * b5 - a1 * b4 + a2 * b3 +
            a3 * b2 - a4 * b1 + a5 * b0;

        if (det == mdp_complex(0)) [[unlikely]]
          error("inv(...)\ndeterminant is zero");

        mdp_complex invdet = mdp_complex(1) / det;

        ans.m_data[0] = (a[5] * b5 - a[6] * b4 + a[7] * b3) * invdet;
        ans.m_data[1] = (-a[1] * b5 + a[2] * b4 - a[3] * b3) * invdet;
        ans.m_data[2] = (a[13] * a5 - a[14] * a4 + a[15] * a3) * invdet;
        ans.m_data[3] = (-a[9] * a5 + a[10] * a4 - a[11] * a3) * invdet;

        ans.m_data[4] = (-a[4] * b5 + a[6] * b2 - a[7] * b1) * invdet;
        ans.m_data[5] = (a[0] * b5 - a[2] * b2 + a[3] * b1) * invdet;
        ans.m_data[6] = (-a[12] * a5 + a[14] * a2 - a[15] * a1) * invdet;
        ans.m_data[7] = (a[8] * a5 - a[10] * a2 + a[11] * a1) * invdet;

        ans.m_data[8] = (a[4] * b4 - a[5] * b2 + a[7] * b0) * invdet;
        ans.m_data[9] = (-a[0] * b4 + a[1] * b2 - a[3] * b0) * invdet;
        ans.m_data[10] = (a[12] * a4 - a[13] * a2 + a[15] * a0) * invdet;
        ans.m_data[11] = (-a[8] * a4 + a[9] * a2 - a[11] * a0) * invdet;

        ans.m_data[12] = (-a[4] * b3 + a[5] * b1 - a[6] * b0) * invdet;
        ans.m_data[13] = (a[0] * b3 - a[1] * b1 + a[2] * b0) * invdet;
        ans.m_data[14] = (-a[12] * a3 + a[13] * a1 - a[14] * a0) * invdet;
        ans.m_data[15] = (a[8] * a3 - a[9] * a1 + a[10] * a0) * invdet;

        return ans;
      }

      // Gauss-Jordan

      mdp_matrix tma(*this);

      std::fill(ans.m_data.begin(), ans.m_data.end(), mdp_complex(0));
      for (mdp_uint i = 0; i < n; i++)
        ans.m_data[i * n + i] = mdp_complex(1);

      for (mdp_uint c = 0; c < n; c++)
      {
        mdp_uint rmax = c;
        mdp_complex pivot = tma(c, c);
        auto pivot_norm = abs2(pivot);

        // choosing pivot

        for (mdp_uint r = c + 1; r < n; r++)
        {
          auto v = tma(r, c);
          auto nv = abs2(v);
          if (nv > pivot_norm)
          {
            rmax = r;
            pivot = v;
            pivot_norm = nv;
          }
        }

        if (pivot == mdp_complex(0)) [[unlikely]]
          error("inv(...)\nsingular matrix");

        std::swap_ranges(&tma(rmax, 0), &tma(rmax, 0) + n, &tma(c, 0));
        std::swap_ranges(&ans(rmax, 0), &ans(rmax, 0) + n, &ans(c, 0));

        mdp_complex invPivot = mdp_complex(1) / tma(c, c);

        for (mdp_uint i = 0; i < n; i++)
        {
          tma(c, i) *= invPivot;
          ans(c, i) *= invPivot;
        }

        // elimination

        for (mdp_uint r = 0; r < n; r++)
        {
          if (r == c)
            continue;

          mdp_complex factor = tma(r, c);

          for (mdp_uint i = c; i < n; i++)
            tma(r, i) -= factor * tma(c, i);

          for (mdp_uint i = 0; i < n; i++)
            ans(r, i) -= factor * ans(c, i);
        }
      }

      return ans;
    }
  };

  inline std::ostream &operator<<(std::ostream &os, const mdp_matrix &a)
  {
    for (mdp_uint i = 0; i < a.rows(); i++)
    {
      if (i == 0)
        os << "[[";
      else
        os << " [";

      for (mdp_uint j = 0; j < a.cols(); j++)
      {
        if (j == 0)
          os << " " << a(i, j);
        else
          os << ", " << a(i, j) << " ";
      }

      if (i == (a.rows() - 1))
        os << "]]\n";
      else
        os << "], \n";
    }
    return os;
  }

  inline bool operator==(const mdp_matrix &a, const mdp_matrix &b)
  {
    if (a.rows() != b.rows() || a.cols() != b.cols())
      return false;

    for (mdp_uint i = 0; i < a.size(); i++)
      if (a[i] != b[i])
        return false;

    return true;
  }

  inline bool operator!=(const mdp_matrix &a, const mdp_matrix &b)
  {
    return !(a == b);
  }

  mdp_matrix operator+(mdp_complex b, const mdp_matrix &a)
  {
    return a + b;
  }

  mdp_matrix operator-(mdp_complex b, const mdp_matrix &a)
  {
    return -(a - b);
  }

  mdp_matrix operator*(mdp_complex x, const mdp_matrix &y)
  {
    return y * x;
  }

  mdp_matrix operator/(mdp_complex b, const mdp_matrix &a)
  {
    return a.inv() * b;
  }

  mdp_matrix operator+(mdp_real b, const mdp_matrix &a)
  {
    return a + b;
  }

  mdp_matrix operator-(mdp_real b, const mdp_matrix &a)
  {
    return -(a - b);
  }

  mdp_matrix operator*(mdp_real a, const mdp_matrix &b)
  {
    return b * a;
  }

  mdp_matrix operator/(mdp_real b, const mdp_matrix &a)
  {
    return a.inv() * b;
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
      {
        tmp(r, c) = (r == c) ? mdp_complex(1) : mdp_complex(0);
      }

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
      {
        tmp(r, c) = mdp_complex(0);
      }

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
    mdp_real current_max = 0.0;

    for (mdp_uint i = 0; i < a.size(); ++i)
      current_max = std::max(current_max, abs(a[i]));

    return current_max;
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
    if (a.rows() < 2 || a.cols() < 2 || i >= a.rows() || j >= a.cols())
      error("submatrix(...)\nWrong dimensions in submatrix");
#endif

    mdp_matrix tmp(a.rows() - 1, a.cols() - 1);

    for (mdp_uint r = 0; r < tmp.rows(); ++r)
    {
      const mdp_uint src_r = (r < i) ? r : r + 1;

      for (mdp_uint c = 0; c < tmp.cols(); ++c)
      {
        const mdp_uint src_c = (c < j) ? c : c + 1;
        tmp(r, c) = a(src_r, src_c);
      }
    }

    return tmp;
  }

  mdp_complex det(const mdp_matrix &a)
  {
    return a.det();
  }

  mdp_matrix inv(const mdp_matrix &a)
  {
    return a.inv();
  }

  mdp_matrix transpose(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.cols(), a.rows());

    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
      {
        tmp(c, r) = a(r, c);
      }

    return tmp;
  }

  mdp_complex trace(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("trace(...)\nmdp_matrix is not squared");
#endif
    mdp_complex x = mdp_complex(0);
    for (mdp_uint c = 0; c < a.cols(); c++)
    {
      x += a(c, c);
    }
    return x;
  }

  mdp_matrix conj(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.rows(), a.cols());
    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
      {
        tmp(r, c) = conj(a(r, c));
      }

    return tmp;
  }

  mdp_matrix hermitian(const mdp_matrix &a)
  {
    mdp_matrix tmp(a.cols(), a.rows());

    // // optimized for 3x3 matrix
    if (a.rows() == 3 && a.cols() == 3)
    {
      const mdp_complex *src = a.address();
      mdp_complex *dst = tmp.address();

      dst[0] = conj(src[0]);
      dst[1] = conj(src[3]);
      dst[2] = conj(src[6]);

      dst[3] = conj(src[1]);
      dst[4] = conj(src[4]);
      dst[5] = conj(src[7]);

      dst[6] = conj(src[2]);
      dst[7] = conj(src[5]);
      dst[8] = conj(src[8]);

      return tmp;
    }

    for (mdp_uint r = 0; r < a.rows(); r++)
      for (mdp_uint c = 0; c < a.cols(); c++)
      {
        tmp(c, r) = conj(a(r, c));
      }

    return tmp;
  }

  mdp_matrix pow(const mdp_matrix &a, int i)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("pow(...)\nmdp_matrix is not squared");
#endif

    const mdp_uint n = a.rows();

    if (i == 0)
      return mdp_identity(n);

    mdp_matrix base = (i < 0) ? inv(a) : a;
    mdp_uint e = (i < 0) ? -i : i;

    mdp_matrix result = mdp_identity(n);

    while (e)
    {
      if (e & 1)
        result = result * base;

      e >>= 1;

      if (e)
        base = base * base;
    }

    return result;
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
      i = -i;
      ++j;
      t1 = ((1.0 * i) / (j)) * c;
      tmp += t1;
    } while (max(t1) > mdp_precision);

    return tmp;
  }

  std::pair<mdp_matrix, mdp_matrix> sincos(const mdp_matrix &a)
  {
#ifdef CHECK_ALL
    if (a.rows() != a.cols())
      error("sincos(...)\nmdp_matrix is not squared");
#endif

    mdp_matrix a2 = a * a;

    mdp_matrix sin_sum = a;
    mdp_matrix cos_sum = mdp_identity(a.rows());

    mdp_matrix sin_t = a;
    mdp_matrix cos_t = mdp_identity(a.rows());

    mdp_uint i_sin = 1;
    mdp_uint i_cos = 0;

    do
    {
      // next sin term
      ++i_sin;
      sin_t = (-1.0 / i_sin) * sin_t * a2;
      ++i_sin;
      sin_t *= 1.0 / i_sin;
      sin_sum += sin_t;

      // next cos term
      ++i_cos;
      cos_t = (-1.0 / i_cos) * cos_t * a2;
      ++i_cos;
      cos_t *= 1.0 / i_cos;
      cos_sum += cos_t;

    } while (max(sin_t) > mdp_precision || max(cos_t) > mdp_precision);

    return {sin_sum, cos_sum};
  }

  mdp_matrix sin(const mdp_matrix &a)
  {
    return sincos(a).first;
  }

  mdp_matrix cos(const mdp_matrix &a)
  {
    return sincos(a).second;
  }
} // namespace MDP

#endif /* MDP_MATRIX_ */
