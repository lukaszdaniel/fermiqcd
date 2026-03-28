#include "mdp.h"

using namespace MDP;

void test_factorial()
{
  assert(mdp_permutations(0) == 1);
  assert(mdp_permutations(1) == 1);
  assert(mdp_permutations(2) == 2);
  assert(mdp_permutations(3) == 6);
  assert(mdp_permutations(4) == 24);
  assert(mdp_permutations(5) == 120);
}

void test_permutation_sort()
{
  {
    int arr[] = {0, 2, 1};
    mdp_permutation_sort(arr, 1);
    assert(arr[0] == 0 && arr[1] == 2 && arr[2] == 1);
  }

  {
    int arr[] = {0, 3, 1};
    mdp_permutation_sort(arr, 2);
    assert(arr[0] == 0 && arr[1] == 1 && arr[2] == 3);
  }

  {
    int arr[] = {2, 1};
    mdp_permutation_sort(arr, 1);
    assert(arr[0] == 1 && arr[1] == 2);
  }
}

void test_permutation_basic()
{
  // Permutations for n = 3:
  // 012 (k = 0)
  // 021 (k = 1)
  // 102 (k = 2)
  // 120 (k = 3)
  // 201 (k = 4)
  // 210 (k = 5)

  std::vector<std::vector<int>> expected = {
      {0, 1, 2},
      {0, 2, 1},
      {1, 0, 2},
      {1, 2, 0},
      {2, 0, 1},
      {2, 1, 0}};

  for (int k = 0; k < 6; k++)
  {
    for (int i = 0; i < 3; i++)
    {
      int val = mdp_permutation(3, k, i);
      assert(val == expected[k][i]);
    }
  }
}

void test_permutation_edges()
{
  // invalid cases
  assert(mdp_permutation(3, 6, 0) == -1); // k == n!
  assert(mdp_permutation(3, 7, 0) == -1); // k > n!
  assert(mdp_permutation(3, 0, 3) == -1); // i == n
  assert(mdp_permutation(3, 0, 4) == -1); // i > n
}

void test_permutation_last()
{
  // last permutation for n = 4 is 3210
  int n = 4;
  int k = mdp_permutations(n) - 1;

  std::vector<int> expected = {3, 2, 1, 0};

  for (int i = 0; i < n; i++)
  {
    assert(mdp_permutation(n, k, i) == expected[i]);
  }
}

void test_permutation_first()
{
  int n = 5;
  int k = 0;

  for (int i = 0; i < n; i++)
  {
    assert(mdp_permutation(n, k, i) == i);
  }
}

int main()
{
  test_factorial();
  test_permutation_sort();
  test_permutation_basic();
  test_permutation_edges();
  test_permutation_first();
  test_permutation_last();

  std::cout << "All tests passed!\n";
  return 0;
}
