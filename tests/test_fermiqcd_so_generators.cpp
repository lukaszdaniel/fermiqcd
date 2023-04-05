#include "fermiqcd.h"

using namespace MDP;

int main()
{
    SO_Generators g(3);
    for (int a = 0; a < g.ngenerators; a++)
        std::cout << "g=" << g.lambda[a] << "\n";
    return 0;
}
