// Program: example08.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  if (mdp.nproc() == 2)
  {
    if (ME == 1)
    {
      int b;
      b = 4 * 8;
      mdp.put(b, 0);
    }
    if (ME == 0)
    {
      int a, b;
      a = 5 * 7;
      mdp.get(b, 1);
      mdp << a + b << "\n";
    }
  }
  else
  {
    mdp << "Sorry, this only runs on 2 processes\n";
  }
  mdp.close_wormholes();
}
