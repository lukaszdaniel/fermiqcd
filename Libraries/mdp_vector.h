/////////////////////////////////////////////////////////////////
/// @file mdp_vector.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_vector
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_VECTOR_
#define MDP_VECTOR_

#include <array>
#include "mdp_global_vars.h"
#define VECTOR_MAX_DIM 10

namespace MDP
{
  /// @brief discrete vectors to navigate on a lattice
  ///
  class mdp_vector
  {
  public:
    using value_type = mdp_uint;
    static constexpr mdp_uint dim = VECTOR_MAX_DIM;

  private:
    std::array<value_type, dim> m_x{};

  public:
    // zero-initialized by default
    constexpr mdp_vector() noexcept = default;

    constexpr mdp_vector(std::initializer_list<value_type> init) noexcept
    {
      std::copy_n(init.begin(),
                  std::min(init.size(), size_t(dim)),
                  m_x.begin());
    }

    // Overloading [] operator to access elements in array style
    constexpr value_type &operator[](std::size_t i) noexcept { return m_x[i]; }
    constexpr value_type operator[](std::size_t i) const noexcept { return m_x[i]; }

    constexpr auto begin() noexcept { return m_x.begin(); }
    constexpr auto end() noexcept { return m_x.end(); }
    constexpr auto begin() const noexcept { return m_x.begin(); }
    constexpr auto end() const noexcept { return m_x.end(); }
  };

  [[nodiscard]]
  constexpr mdp_vector binary2versor(mdp_uint a) noexcept
  {
    mdp_vector v;
    for (std::size_t i = 0; i < mdp_vector::dim; ++i)
      v[i] = (a >> i) & 0x1;
    return v;
  }

  [[nodiscard]]
  constexpr mdp_uint vector2binary(const mdp_vector &v)
  {
#ifdef CHECK_ALL
    for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    {
      if (v[i] != 0 && v[i] != 1)
        error("vector2binary");
    }
#endif

    mdp_uint result = 0;
    for (std::size_t i = 0; i < mdp_vector::dim; ++i)
      result |= (static_cast<mdp_uint>(v[i]) << i);

    return result;
  }
} // namespace MDP

#endif /* MDP_VECTOR_ */
