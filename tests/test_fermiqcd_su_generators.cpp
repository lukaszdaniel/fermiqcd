#include "fermiqcd.h"

int main()
{
    SU_Generators g(3);
    for (int a = 0; a < g.ngenerators; a++)
        std::cout << "g=" << g.lambda[a] << "\n";
    return 0;
}