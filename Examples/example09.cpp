// Program: example09.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int a;
  if (mdp.nproc() == 2)
  {
    if (isSubProcess(1))
    {
      a = 4 * 8;
      mdp.add(a);
    }
    if (isMainProcess())
    {
      a = 5 * 7;
      mdp.add(a);
      std::cout << "a=" << a << "\n";
    }
  }
  else
  {
    mdp << "Sorry, this only runs on 2 processes\n";
  }
  mdp.close_wormholes();
}
