#include "mdp.h"

using namespace MDP;

int main()
{
    mdp_measure m;
    // store 10 measurements
    for (int i = 0; i < 10; i++)
        m << 3.0 + mdp_random.gaussian(2.0);

    std::cout << m.getmean() << " +/- " << m.getmerr() << std::endl;
    std::cout << m.getnum() << std::endl;

    m = sin(exp(m) + m);

    std::cout << m.getmean() << " +/- " << m.getmerr() << std::endl;
    std::cout << m.getnum() << std::endl;

    return 0;
}
