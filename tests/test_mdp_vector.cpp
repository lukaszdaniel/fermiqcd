#include "mdp.h"

using namespace MDP;

void test_default_constructor()
{
  mdp_vector v;

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    assert(v[i] == 0);
}

void test_initializer_list()
{
  mdp_vector v{1, 2, 3};

  assert(v[0] == 1);
  assert(v[1] == 2);
  assert(v[2] == 3);

  // reszta powinna być 0
  for (std::size_t i = 3; i < mdp_vector::dim; ++i)
    assert(v[i] == 0);
}

void test_initializer_truncation()
{
  // więcej elementów niż dim
  mdp_vector v{
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
      11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    assert(v[i] == static_cast<mdp_uint>(i + 1));
}

void test_operator_access()
{
  mdp_vector v;

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    v[i] = static_cast<mdp_uint>(i * 2);

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    assert(v[i] == static_cast<mdp_uint>(i * 2));
}

void test_const_access()
{
  const mdp_vector v{5, 6, 7};

  assert(v[0] == 5);
  assert(v[1] == 6);
  assert(v[2] == 7);
}

void test_iterators()
{
  mdp_vector v;
  std::iota(v.begin(), v.end(), 0);

  mdp_uint sum = std::accumulate(v.begin(), v.end(), 0);

  mdp_uint expected = (mdp_vector::dim - 1) * mdp_vector::dim / 2;
  assert(sum == expected);
}

// ===================== binary2versor =====================

void test_binary2versor_basic()
{
  mdp_vector v = binary2versor(0b101);

  assert(v[0] == 1);
  assert(v[1] == 0);
  assert(v[2] == 1);
}

void test_binary2versor_zero()
{
  mdp_vector v = binary2versor(0);

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    assert(v[i] == 0);
}

void test_binary2versor_full()
{
  mdp_uint all_ones = 0;
  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    all_ones |= (static_cast<mdp_uint>(1) << i);

  mdp_vector v = binary2versor(all_ones);

  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    assert(v[i] == 1);
}

// ===================== vector2binary =====================

void test_vector2binary_basic()
{
  mdp_vector v{1, 0, 1};

  mdp_uint result = vector2binary(v);

  assert(result == 0b101);
}

void test_vector2binary_zero()
{
  mdp_vector v;

  assert(vector2binary(v) == 0);
}

void test_vector2binary_full()
{
  mdp_vector v;
  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    v[i] = 1;

  mdp_uint expected = 0;
  for (std::size_t i = 0; i < mdp_vector::dim; ++i)
    expected |= (static_cast<mdp_uint>(1) << i);

  assert(vector2binary(v) == expected);
}

// ===================== round-trip =====================

void test_round_trip()
{
  // test dla kilku wartości
  for (mdp_int x = 0; x < (1 << std::min<std::size_t>(mdp_vector::dim, 10)); ++x)
  {
    mdp_vector v = binary2versor(x);
    mdp_int y = vector2binary(v);

    assert(x == y);
  }
}

// ===================== constexpr tests =====================

constexpr bool test_constexpr()
{
  constexpr mdp_vector v = binary2versor(0b1011);
  static_assert(v[0] == 1);
  static_assert(v[1] == 1);
  static_assert(v[2] == 0);
  static_assert(v[3] == 1);

  constexpr mdp_uint x = vector2binary(v);
  static_assert(x == 0b1011);

  return true;
}

static_assert(test_constexpr());

// ===================== MAIN =====================

int main()
{
  test_default_constructor();
  test_initializer_list();
  test_initializer_truncation();
  test_operator_access();
  test_const_access();
  test_iterators();

  test_binary2versor_basic();
  test_binary2versor_zero();
  test_binary2versor_full();

  test_vector2binary_basic();
  test_vector2binary_zero();
  test_vector2binary_full();

  test_round_trip();

  std::cout << "All mdp_vector tests passed!\n";
}
