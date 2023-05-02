#include <algorithm>
#include "mdp.h"
#include "MDP_PVegas.h"

using namespace MDP;

class MyFunction : public VegasBase
{
protected:
  /** @brief Sample function to be integrated
   *
   * @note its exact integral formula is:
   * I =  c1 * y + c2 - cos(x) sin(y)
   *
   * @note exact value is:
   * val = std::sin(1) * (1 - std::cos(1))
   */
  double m_function(const double *x)
  {
    return std::sin(x[0]) * std::cos(x[1]);
  }
};

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  MyFunction myfunction;
  mdp << "Exact value = " << std::sin(1) * (1 - std::cos(1)) << "\n";
  myfunction.setDimension(2);
  myfunction.setIntegrationLimits(0, 0, 1);
  myfunction.setIntegrationLimits(1, 0, 1);
  mdp << "Intergal = " << myfunction.Integrate();
  mdp.close_wormholes();
  return 0;
}
