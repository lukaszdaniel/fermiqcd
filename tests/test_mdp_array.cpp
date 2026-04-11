#include <iostream>
#include <cassert>
#include "mdp.h"

using namespace MDP;

// --- Funkcje pomocnicze do applytoall ---
float add_const(float x)
{
  return x + 5;
}

float add_two(float a, float b)
{
  return a + b;
}

// --- TESTY ---
int main()
{
  std::cout << "==== TEST mdp_array ====\n";

  // 1. Konstruktor + operator()
  mdp_array<float, 3> a(3, 3, 3);
  for (mdp_uint i = 0; i < 3; i++)
    for (mdp_uint j = 0; j < 3; j++)
      for (mdp_uint k = 0; k < 3; k++)
        a(i, j, k) = i * 100 + j * 10 + k;

  std::cout << "Array a:\n"
            << a << std::endl;

  // 2. operator[]
  std::cout << "Test operator[]:\n";
  for (mdp_uint i = 0; i < a.size(); i++)
    std::cout << a[i] << " ";
  std::cout << "\n";

  // 3. Kopiowanie (copy constructor)
  mdp_array<float, 3> b = a;
  std::cout << "Array b (copy of a):\n"
            << b << std::endl;

  // 4. Przypisanie
  mdp_array<float, 3> c(3, 3, 3);
  c = a;
  std::cout << "Array c (assigned from a):\n"
            << c << std::endl;

  // 5. operator +
  mdp_array<float, 3> d = a + b;
  std::cout << "a + b:\n"
            << d << std::endl;

  // 6. operator -
  mdp_array<float, 3> e = d - a;
  std::cout << "d - a:\n"
            << e << std::endl;

  // 7. operator * (skalar)
  mdp_array<float, 3> f = 2.0f * a;
  std::cout << "2 * a:\n"
            << f << std::endl;

  // 8. applytoall (unarny)
  mdp_array<float, 3> g = applytoall(a, add_const);
  std::cout << "a + 5 (applytoall):\n"
            << g << std::endl;

  // 9. applytoall (binarny)
  mdp_array<float, 3> h = applytoall(a, b, add_two);
  std::cout << "a + b (applytoall binary):\n"
            << h << std::endl;

  // 10. Test dims() i size()
  const auto &dims = a.dims();
  std::cout << "Dims a: ";
  for (mdp_suint i = 0; i < a.ndim(); i++)
    std::cout << dims[i] << " ";
  std::cout << "\nTotal size: " << a.size() << "\n";

  // 11. Test address()
  float *ptr = a.address();
  std::cout << "First element via pointer: " << ptr[0] << "\n";

  // 12. Test dimension() (realloc)
  mdp_uint newdims[3] = {2, 2, 2};
  a.dimension(newdims);
  std::cout << "After resizing a to 2x2x2, size = " << a.size() << "\n";

  // 13. Prosty assert
  mdp_array<int, 1> test(5);
  for (int i = 0; i < 5; i++)
    test[i] = i;

  for (int i = 0; i < 5; i++)
    assert(test[i] == i);

  // Other dimensions:
  mdp_array<float, 1> a1(3);
  for (mdp_uint i = 0; i < 3; i++)
    a1(i) = i;

  std::cout << "Array a1:\n"
            << a1 << std::endl;

  mdp_array<float, 2> a2(3, 3);
  for (mdp_uint i = 0; i < 3; i++)
    for (mdp_uint j = 0; j < 3; j++)
      a2(i, j) = i * 10 + j;

  std::cout << "Array a2:\n"
            << a2 << std::endl;

  std::cout << "Assert test passed!\n";

  std::cout << "==== ALL TESTS DONE ====\n";
  return 0;
}
