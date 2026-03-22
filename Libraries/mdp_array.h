/////////////////////////////////////////////////////////////////
/// @file mdp_array.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_array
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_ARRAY_
#define MDP_ARRAY_

#include <iostream>
#include <memory>
#include <array>
#include <numeric>
#include <span>
#include "mdp_macros.h"
#include "mdp_global_vars.h"
#define ARRAY_MAX_DIM 5

namespace MDP
{
  /// @brief generic container for multidimensional arrays
  ///
  /// Example:
  /// @verbatim
  ///    mdp_array<mdp_real,3> a(5,5,5);
  ///    a(0,0,0)=3.15;
  /// @endverbatim
  template <class T, mdp_uint ndims>
  class mdp_array
  {
  private:
    using Array = std::array<mdp_uint, ARRAY_MAX_DIM>;
    Array m_size{}; // size of each dimension

    std::unique_ptr<T[]> m_owner;
    std::span<T> m_data;

    /** @brief Assign dimension values to this array
     *
     * @param c0_ size of 0-th dimension
     * @param c1_ size of 1-st dimension
     * @param c2_ size of 2-nd dimension
     * @param c3_ size of 3-rd dimension
     * @param c4_ size of 4-th dimension
     */
    void calculateDimensions(const mdp_uint c0_, const mdp_uint c1_,
                             const mdp_uint c2_, const mdp_uint c3_, const mdp_uint c4_)
    {
      std::array<mdp_uint, 5> args = {c0_, c1_, c2_, c3_, c4_};
      for (mdp_uint i = 0; i < ARRAY_MAX_DIM; i++)
        m_size[i] = (i < ndims) ? args[i] : 1;
    }

    /** @brief Assign dimension values to this array
     *
     * @param p Reference to the array describing dimensions
     */
    void calculateDimensions(const Array &p)
    {
      for (mdp_uint i = 0; i < ndims; i++)
        m_size[i] = p[i];
      for (mdp_uint i = ndims; i < ARRAY_MAX_DIM; i++)
        m_size[i] = 1;
    }

    /** @brief Assign dimension values to this array from raw pointer
     *
     * @param p Pointer to array of dimension sizes
     */
    void calculateDimensions(const mdp_uint *p)
    {
      for (mdp_uint i = 0; i < ndims; i++)
        m_size[i] = p[i];
      for (mdp_uint i = ndims; i < ARRAY_MAX_DIM; i++)
        m_size[i] = 1;
    }

    /** @brief Allocate required memory
     */
    void allocate()
    {
      auto n = total_size();
      if (n == 0)
      {
        m_owner.reset();
        m_data = {};
        return;
      }

      m_owner = std::make_unique<T[]>(n);
      m_data = std::span<T>(m_owner.get(), n);
    }

    void bind_external(T *ptr)
    {
      m_owner.reset();
      m_data = std::span<T>(ptr, total_size());
    }

    void validateIndices(const mdp_uint i0, const mdp_uint i1, const mdp_uint i2,
                         const mdp_uint i3, const mdp_uint i4) const
    {
      if (i0 >= m_size[0])
        error("mdp_array::operator()\nIndex out of bounds (dimension 0)");
      if constexpr (ndims >= 2)
      {
        if (i1 >= m_size[1])
          error("mdp_array::operator()\nIndex out of bounds (dimension 1)");
      }
      if constexpr (ndims >= 3)
      {
        if (i2 >= m_size[2])
          error("mdp_array::operator()\nIndex out of bounds (dimension 2)");
      }
      if constexpr (ndims >= 4)
      {
        if (i3 >= m_size[3])
          error("mdp_array::operator()\nIndex out of bounds (dimension 3)");
      }
      if constexpr (ndims >= 5)
      {
        if (i4 >= m_size[4])
          error("mdp_array::operator()\nIndex out of bounds (dimension 4)");
      }
    }

    /** Copy one array to another array
     *
     * @param a Array to be copied
     */
    void copyArray(const mdp_array &a)
    {
      m_size = a.m_size;
      allocate();
      std::copy(a.m_data.begin(), a.m_data.end(), m_data.begin());
    }

    mdp_uint total_size() const
    {
      return std::accumulate(m_size.begin(), m_size.end(), mdp_uint{1}, std::multiplies<>());
    }

  public:
    /** @brief Dimension of this array
     *
     * @return dimension of this array
     *
     * @note The value is at most equal to 5
     */
    [[nodiscard]] mdp_uint ndim() const { return ndims; }

    /** @brief Begining of this array
     *
     * @return Pointer to the begining of this array
     */
    [[nodiscard]] T *address() { return m_data.data(); }
    [[nodiscard]] const T *address() const { return m_data.data(); }

    // C++20: iterators for range-based for loops
    [[nodiscard]] auto begin() noexcept { return m_data.begin(); }
    [[nodiscard]] auto end() noexcept { return m_data.end(); }
    [[nodiscard]] auto begin() const noexcept { return m_data.begin(); }
    [[nodiscard]] auto end() const noexcept { return m_data.end(); }

    /** @brief Size of each dimension
     *
     * @return Pointer to the array of sizes
     */
    [[nodiscard]] const Array &dims() const { return m_size; }

    /** @brief Size of i-th dimension
     *
     * @return Size of i-th dimension
     */
    [[nodiscard]] mdp_uint size(mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= ndims)
        error("mdp_array::size(...)\nIndex out of bounds");
#endif
      return m_size[i];
    }

    /** @brief Total size of this array
     *
     * @return Total number of elements
     *
     * @note This is a product of all dimensions
     */
    [[nodiscard]] mdp_uint size() const { return static_cast<mdp_uint>(m_data.size()); }

    T &operator[](mdp_uint i)
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_array::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    const T &operator[](mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_array::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    void dimension(const mdp_uint *p)
    {
      calculateDimensions(p);
      allocate();
    }

    void dimension(const mdp_uint c0_ = 1, const mdp_uint c1_ = 1,
                   const mdp_uint c2_ = 1, const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
      allocate();
    }

    mdp_array(const mdp_uint c0_ = 1, const mdp_uint c1_ = 1,
              const mdp_uint c2_ = 1, const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
      allocate();
    }

    mdp_array(const Array &dims)
    {
      calculateDimensions(dims);
      allocate();
    }

    mdp_array(T *ptr, const mdp_uint c0_ = 1, const mdp_uint c1_ = 1,
              const mdp_uint c2_ = 1, const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
      bind_external(ptr);
    }

    mdp_array(const mdp_array &a)
    {
      copyArray(a);
    }

    ~mdp_array() = default;
    mdp_array(mdp_array &&) noexcept = default;
    mdp_array &operator=(mdp_array &&) noexcept = default;

    mdp_array &operator=(const mdp_array &a)
    {
      if (this != &a)
        copyArray(a);
      return *this;
    }

    T &operator()(mdp_uint i0, mdp_uint i1 = 0, mdp_uint i2 = 0,
                  mdp_uint i3 = 0, mdp_uint i4 = 0)
    {
#ifdef CHECK_ALL
      validateIndices(i0, i1, i2, i3, i4);
#endif
      return m_data[(((i0 * m_size[1] + i1) * m_size[2] + i2) * m_size[3] + i3) * m_size[4] + i4];
    }

    const T &operator()(mdp_uint i0, mdp_uint i1 = 0, mdp_uint i2 = 0,
                        mdp_uint i3 = 0, mdp_uint i4 = 0) const
    {
#ifdef CHECK_ALL
      validateIndices(i0, i1, i2, i3, i4);
#endif
      return m_data[(((i0 * m_size[1] + i1) * m_size[2] + i2) * m_size[3] + i3) * m_size[4] + i4];
    }
  };

  template <class T, mdp_uint ndims>
  mdp_array<T, ndims> operator+(const mdp_array<T, ndims> &a, const mdp_array<T, ndims> &b)
  {
    if (a.dims() != b.dims())
      error("operator+: incompatible dimensions");
    mdp_array<T, ndims> tmp(a.dims());
    std::transform(a.begin(), a.end(),
                   b.begin(), tmp.begin(), std::plus<T>{});
    return tmp;
  }

  template <class T, mdp_uint ndims>
  mdp_array<T, ndims> operator-(const mdp_array<T, ndims> &a, const mdp_array<T, ndims> &b)
  {
    if (a.dims() != b.dims())
      error("operator-: incompatible dimensions");
    mdp_array<T, ndims> tmp(a.dims());
    std::transform(a.begin(), a.end(),
                   b.begin(), tmp.begin(), std::minus<T>{});
    return tmp;
  }

  template <class T, class T2, mdp_uint ndims>
  mdp_array<T, ndims> operator*(T2 x, const mdp_array<T, ndims> &a)
  {
    mdp_array<T, ndims> tmp(a.dims());
    std::transform(a.begin(), a.end(),
                   tmp.begin(), [x](T val)
                   { return val * x; });
    return tmp;
  }

  // Implementation of generic unary operator
  template <class T, mdp_uint ndims, class Func>
  mdp_array<T, ndims> applytoall(const mdp_array<T, ndims> &a, Func &&f)
  {
    mdp_array<T, ndims> tmp(a.dims());
    std::transform(a.begin(), a.end(),
                   tmp.begin(), std::forward<Func>(f));
    return tmp;
  }

  template <class T, mdp_uint ndims, class Func>
  mdp_array<T, ndims> applytoall(const mdp_array<T, ndims> &a, const mdp_array<T, ndims> &b, Func &&f)
  {
    if (a.dims() != b.dims())
      error("applytoall: incompatible dimensions");
    mdp_array<T, ndims> tmp(a.dims());
    std::transform(a.begin(), a.end(),
                   b.begin(), tmp.begin(), std::forward<Func>(f));
    return tmp;
  }

  // Helper function for pretty printing
  template <class T, mdp_uint ndims>
  void print_recursive_pretty(std::ostream &os,
                              const mdp_array<T, ndims> &a,
                              mdp_uint dim,
                              mdp_uint offset,
                              const std::array<mdp_uint, ARRAY_MAX_DIM> &strides,
                              mdp_uint indent)
  {
    const auto &dims = a.dims();
    mdp_uint current_dim_size = dims[dim];

    std::string indent_str(indent, ' ');
    std::string next_indent_str(indent + 2, ' ');

    // print most inner dimension in one line
    if (dim == ndims - 1)
    {
      os << "{ ";
      for (mdp_uint i = 0; i < current_dim_size; i++)
      {
        os << a[offset + i * strides[dim]];
        if (i != current_dim_size - 1)
          os << ", ";
      }
      os << " }";
      return;
    }

    // other dimensions in block
    os << "{\n";

    for (mdp_uint i = 0; i < current_dim_size; i++)
    {
      os << next_indent_str;

      print_recursive_pretty(os,
                             a,
                             dim + 1,
                             offset + i * strides[dim],
                             strides,
                             indent + 2);

      if (i != current_dim_size - 1)
        os << ",";

      os << "\n";
    }

    os << indent_str << "}";
  }

  template <class T, mdp_uint ndims>
  std::ostream &operator<<(std::ostream &os, const mdp_array<T, ndims> &a)
  {
    if (ndims == 0 || a.size() == 0)
    {
      os << "{}";
      return os;
    }

    std::array<mdp_uint, ARRAY_MAX_DIM> strides;
    strides.fill(1);
    const auto &dims = a.dims();

    if constexpr (ndims > 1)
    {
      strides[ndims - 1] = 1;
      for (int i = static_cast<int>(ndims) - 2; i >= 0; i--)
        strides[i] = strides[i + 1] * dims[i + 1];
    }

    print_recursive_pretty(os, a, 0, 0, strides, 0);

    return os;
  }
} // namespace MDP

#endif /* MDP_ARRAY_ */
