#include "mdp.h"

using namespace MDP;

int main()
{
    mdp_postscript ps("test.eps");
    ps.size(0, 0, 595, 842);
    ps.color(0.2, 0.2, 0.7);
    ps.line(171, 400, 441, 450);
    ps.font("Times-Roman", 12);
    ps.print(441, 450, "a line from (0,0) to here");
    ps.box(120, 350, 170, 394, 1);
    ps.circle(250, 425, 20, 1);

    return 0;
}