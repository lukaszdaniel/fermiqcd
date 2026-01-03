#include "fermiqcd.h"

using namespace MDP;

// Tests equivalence between old path functions and new Path class

int main()
{
  constexpr int size1 = 2;
  constexpr int size2 = 4;

  constexpr int length = 2 * size1 + 2 * size2;
  int old_path[length][2];
  Path path(length);
  Path path2({{1, 2}, {1, 2}, {1, 3}, {1, 3}, {1, 3}, {1, 3}, {-1, 2}, {-1, 2}, {-1, 3}, {-1, 3}, {-1, 3}, {-1, 3}});

  // make a generic path
  for (int i = 0; i < size1; i++)
  {
    path[i][0] = +1;
    path[i + size1 + size2][0] = -1;
    old_path[i][0] = +1;
    old_path[i + size1 + size2][0] = -1;
  }

  for (int i = size1; i < size1 + size2; i++)
  {
    path[i][0] = +1;
    path[i + size1 + size2][0] = -1;
    old_path[i][0] = +1;
    old_path[i + size1 + size2][0] = -1;
  }

  // loop over all possible paths
  for (int mu = 1; mu < 4; mu++)
  {
    for (int nu = mu + 1; nu < 4; nu++)
    {
      // build each path
      for (int i = 0; i < size1; i++)
      {
        path[i][1] = path[i + size1 + size2][1] = mu;
        old_path[i][1] = old_path[i + size1 + size2][1] = mu;
      }
      for (int i = size1; i < size1 + size2; i++)
      {
        path[i][1] = path[i + size1 + size2][1] = nu;
        old_path[i][1] = old_path[i + size1 + size2][1] = nu;
      }
    }
  }

  std::cout << "old_path, path, path2" << std::endl;
  for (int i = 0; i < length; i++)
  {
    std::cout << i << " {" << old_path[i][0] << ", " << old_path[i][1] << "}";
    std::cout << " {" << path[i][0] << ", " << path[i][1] << "}";
    std::cout << " {" << path2[i][0] << ", " << path2[i][1] << "} " << std::endl;
  }

  std::cout << std::endl;

  constexpr int mu = 0;
  constexpr int nu = 1;
  const Path base_path({{1, mu}, {1, mu}, {1, nu}, {-1, mu}, {-1, mu}, {-1, nu}});
  std::cout << "path      0 = ";
  for (auto &[orientation, direction] : base_path)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;

  Path path_90 = base_path;
  path_90.rotate_path(90, mu, nu);
  std::cout << "path     90 = ";
  for (auto &[orientation, direction] : path_90)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  Path path_180 = base_path;
  path_180.rotate_path(180, mu, nu);
  path_90.rotate_path(90, mu, nu);
  std::cout << "path 2 x 90 = ";
  for (auto &[orientation, direction] : path_90)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << "path    180 = ";
  for (auto &[orientation, direction] : path_180)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  Path path_270 = base_path;
  path_270.rotate_path(270, mu, nu);
  path_180.rotate_path(90, mu, nu);
  std::cout << "path 3 x 90 = ";
  for (auto &[orientation, direction] : path_180)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << "path    270 = ";
  for (auto &[orientation, direction] : path_270)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  Path path_360 = base_path;
  path_360.rotate_path(270, mu, nu);
  path_360.rotate_path(90, mu, nu);
  path_270.rotate_path(90, mu, nu);
  std::cout << "path 4 x 90 = ";
  for (auto &[orientation, direction] : path_270)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  std::cout << "path    360 = ";
  for (auto &[orientation, direction] : path_360)
  {
    std::cout << "{" << orientation << ", " << (direction == 0 ? "mu" : "nu") << "}, ";
  }
  std::cout << std::endl;
  return 0;
}
