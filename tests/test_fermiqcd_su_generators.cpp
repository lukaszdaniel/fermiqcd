#include "fermiqcd.h"

using namespace MDP;

int main()
{
    SU_Generators g(3);
    for (mdp_uint a = 0; a < g.ngenerators; a++)
        std::cout << "g(" << a << "):\n" << g.lambda(a) << "\n";
    return 0;
}
