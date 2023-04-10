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

#include <ostream>
#include <memory>
#include <array>
#include <numeric>
#define ARRAY_MAX_DIM 5

namespace MDP
{
  /// @brief generic container for multidimensional arrays
  ///
  /// Example:
  /// @verbatim
  ///    mdp_array<float,3> a(5,5,5);
  ///    a(0,0,0)=3.15;
  /// @endverbatim
  template <class T, mdp_uint ndims>
  class mdp_array
  {
  private:
    bool m_shared;
    mdp_uint m_dim;                             // dimension of the array
    std::array<mdp_uint, ARRAY_MAX_DIM> m_size; // size of each dimension
    std::unique_ptr<T[]> m_data;                // data

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
      m_size = {c0_, c1_, c2_, c3_, c4_};
    }

    /** @brief Assign dimension values to this array
     *
     * @param p Pointer to array describing dimensions
     */
    void calculateDimensions(const mdp_uint *p)
    {
      for (mdp_uint i = 0; i < m_dim; i++)
        m_size[i] = p[i];
      for (mdp_uint i = m_dim; i < ARRAY_MAX_DIM; i++)
        m_size[i] = 1;
    }

    /** @brief Allocate required memory
     */
    void allocate()
    {
      if (!m_shared)
      {
        m_data = std::make_unique<T[]>(size());
        if (size() && !m_data)
          error("mdp_array::allocate()\nOut of memory");
        // memset(m_data, 0, size() * sizeof(T));
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

    void validateIndices(const mdp_uint i0, const mdp_uint i1, const mdp_uint i2,
                         const mdp_uint i3, const mdp_uint i4)
    {
      if ((i1 != 0 && m_dim < 2) ||
          (i2 != 0 && m_dim < 3) ||
          (i3 != 0 && m_dim < 4) ||
          (i4 != 0 && m_dim < ARRAY_MAX_DIM))
        error("mdp_array::operator()(...)\nIncompatible size()");
      if (i0 >= m_size[0])
        error("mdp_array::operator()\nIndex out of bounds");
      if (i1 >= m_size[1])
        error("mdp_array::operator()\nIndex out of bounds");
      if (i2 >= m_size[2])
        error("mdp_array::operator()\nIndex out of bounds");
      if (i3 >= m_size[3])
        error("mdp_array::operator()\nIndex out of bounds");
      if (i4 >= m_size[4])
        error("mdp_array::operator()\nIndex out of bounds");
    }

    /** Copy one array to another array
     *
     * @param a Array to be copied
     */
    void copyArray(const mdp_array &a)
    {
      // if (size() != a.size())
      dimension(a.size());

      for (mdp_uint i = 0; i < ARRAY_MAX_DIM; i++)
        m_size[i] = a.m_size[i];
      for (mdp_uint i = 0; i < size(); i++)
        m_data[i] = a.m_data[i];
    }

  public:
    /** @brief Dimension of this array
     *
     * @return dimension of this array
     *
     * @note The value is at most equal to 5
     */
    mdp_uint ndim() const
    {
      return m_dim;
    }

    /** @brief Begining of this array
     *
     * @return Pointer to the begining of this array
     */
    T *address()
    {
      return m_data.get();
    }

    /** @brief Size of each dimension
     *
     * @return Pointer to the array of sizes
     */
    const mdp_uint *dims() const
    {
      return m_size.data();
    }

    T &operator[](const mdp_uint i)
    {
#ifdef CHECK_ALL
      if (i >= size())
        error("mdp_array::operator[]\nIndex out of bounds");
#endif
      return m_data[i];
    }

    const T &operator[](const mdp_uint i) const
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
    mdp_uint size(mdp_uint i) const
    {
#ifdef CHECK_ALL
      if (i >= size())
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
    mdp_uint size() const
    {
      return std::accumulate(m_size.begin(), m_size.end(), 1, std::multiplies<mdp_uint>());
    }

    void dimension(const mdp_uint *p)
    {
      calculateDimensions(p);
      reallocate();
    }

    void dimension(const mdp_uint c0_ = 1, const mdp_uint c1_ = 1,
                   const mdp_uint c2_ = 1, const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
      allocate();
    }

    mdp_array(const mdp_uint c0_ = 1, const mdp_uint c1_ = 1, const mdp_uint c2_ = 1,
              const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
        : m_shared(false), m_dim(ndims), m_data(nullptr)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
      allocate();
    }

    mdp_array(const mdp_uint *p)
        : m_shared(false), m_dim(ndims), m_data(nullptr)
    {
      calculateDimensions(p);
      allocate();
    }

    mdp_array(const T *m0, const mdp_uint c0_ = 1, const mdp_uint c1_ = 1,
              const mdp_uint c2_ = 1, const mdp_uint c3_ = 1, const mdp_uint c4_ = 1)
        : m_shared(true), m_data(m0), m_dim(ndims)
    {
      calculateDimensions(c0_, c1_, c2_, c3_, c4_);
    }

    mdp_array(const T *m0, const mdp_uint *p)
        : m_shared(true), m_data(m0), m_dim(ndims)
    {
      calculateDimensions(p);
    }

    mdp_array(const mdp_array &a)
        : m_shared(false), m_dim(ndims), m_data(nullptr)
    {
      if (m_dim != a.m_dim)
        error("mdp_array::mdp_array(...)\nIncompatible size()");

      copyArray(a);
    }

    virtual ~mdp_array()
    {
      deallocate();
    }

    const mdp_array &operator=(const mdp_array &a)
    {
      m_dim = a.m_dim;
      copyArray(a);
      return *this;
    }

    T &operator()(const mdp_uint i0, const mdp_uint i1 = 0, const mdp_uint i2 = 0,
                  const mdp_uint i3 = 0, const mdp_uint i4 = 0)
    {
#ifdef CHECK_ALL
      validateIndices(i0, i1, i2, i3, i4);
#endif
      return m_data[(((i0 * m_size[1] + i1) * m_size[2] + i2) * m_size[3] + i3) * m_size[4] + i4];
    }

    const T &operator()(const mdp_uint i0, const mdp_uint i1 = 0, const mdp_uint i2 = 0,
                        const mdp_uint i3 = 0, const mdp_uint i4 = 0) const
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
    mdp_array<T, ndims> tmp(a.dims());
    for (mdp_uint i = 0; i < a.size(); i++)
      tmp[i] = a[i] + b[i];
    return tmp;
  }

  template <class T, mdp_uint ndims>
  mdp_array<T, ndims> operator-(const mdp_array<T, ndims> &a, const mdp_array<T, ndims> &b)
  {
    mdp_array<T, ndims> tmp(a.dims());
    for (mdp_uint i = 0; i < a.size(); i++)
      tmp[i] = a[i] - b[i];
    return tmp;
  }

  template <class T, class T2, mdp_uint ndims>
  mdp_array<T, ndims> operator*(T2 x, const mdp_array<T, ndims> &a)
  {
    mdp_array tmp(a.dims());
    for (mdp_uint i = 0; i < a.size(); i++)
      tmp[i] = a[i] * x;
    return tmp;
  }

  // Implementation of generic unary operator
  template <class T, mdp_uint ndims>
  mdp_array<T, ndims> applytoall(const mdp_array<T, ndims> &a, T (*fptr)(T, void *),
                                 void *x = nullptr)
  {
    mdp_array<T, ndims> tmp(a.dims());
    for (mdp_uint i = 0; i < a.size(); i++)
      tmp[i] = (*fptr)(a[i], x);
    return tmp;
  }

  // implementation of generic binary operator
  template <class T, mdp_uint ndims>
  mdp_array<T, ndims> applytoall(const mdp_array<T, ndims> &a, const mdp_array<T, ndims> &b,
                                 T (*fptr)(T, T, void *), void *x = nullptr)
  {
    mdp_array tmp(a.dims());
    for (mdp_uint i = 0; i < a.size(); i++)
      tmp[i] = (*fptr)(a[i], b[i], x);
    return tmp;
  }

  template <class T, mdp_uint ndims>
  std::ostream &operator<<(std::ostream &os, const mdp_array<T, ndims> &a)
  {
    switch (a.ndim())
    {
    case 0:
      os << "{}";
      break;
    case 1:
      os << "{";
      for (mdp_uint i = 0; i < a.length(); i++)
        if (i == 0)
          os << " " << a[i];
        else
          os << ", " << a[i];
      os << "}";
      break;
    case 2:
      for (mdp_uint i = 0; i < a.size(0); i++)
      {
        if (i == 0)
          os << "{{";
        else
          os << " {";
        for (mdp_uint j = 0; j < a.size(1); j++)
          if (j == 0)
            os << " " << a(i, j);
          else
            os << ", " << a(i, j) << " ";
        if (i == (a.size(0) - 1))
          os << "}}\n";
        else
          os << "}, \n";
      }
      break;
    default:
      os << " Sorry. Option not implemented";
      break;
    }

    return os;
  }
} // namespace MDP

#endif /* MDP_ARRAY_ */
