// Program: example13.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int j = 0;
  for (int i = 0; i < 5; i++)
  {
    if (isMainProcess())
      j = i;
    else
      j = 0;

    if (i % 2 == 0)
      mdp.barrier();

    mdp.broadcast(j, 0);
    std::cout << "I am process " << ME
         << ", i=" << i
         << ", j=" << j << "\n";
  }
  mdp.close_wormholes();
}
