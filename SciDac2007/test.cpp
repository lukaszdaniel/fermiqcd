/*
    python3 iqcd.py code +hot TxXxYxZ=4x4x4x4 +loop n=10 { +heatbath +plaquette } +loop n=10 { +ape_smear steps=1 +plaquette +topological_charge }

    +hot warning: assuming default argument nc=3
    +heatbath warning: assuming default argument beta=5.0
    +heatbath warning: assuming default argument steps=1
    +ape_smear warning: assuming default argument alpha=0.7
    +ape_smear warning: assuming default argument cooling_steps=10
    +topological_charge warning: assuming default argument filename=topological_charge_*.vtk
    +topological_charge warning: assuming default argument t=-1
*/
#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
   mdp.open_wormholes(argc, argv);
   std::string filename;
   coefficients coeff;
   constexpr Box L = {4, 4, 4, 4};
   mdp_lattice spacetime(L);
   int nc = 3;
   gauge_field U(spacetime, nc);
   set_hot(U);

   for (int i0 = 0; i0 < 10; i0++)
   {
      coeff["beta"] = 5.0;
      WilsonGaugeAction::heatbath(U, coeff, 1);
      mdp << "average plaquette = " << average_plaquette(U) << "\n";
   }

   for (int i0 = 0; i0 < 10; i0++)
   {
      ApeSmearing::smear(U, 0.7, 1, 10);
      mdp << "average plaquette = " << average_plaquette(U) << "\n";
      {
         float tc = topological_charge_vtk(U, "topological_charge_*.vtk", -1);
         mdp << "total topological charge = " << tc << "\n";
      }
   }

   mdp.close_wormholes();
   return 0;
}
