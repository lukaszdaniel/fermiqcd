#include "mdp.h"

using namespace MDP;

// Example of a program that uses the Leverberger-Marquardt

// fitting function
float f(float x, float *a, mdp_int ma, void *junk)
{
    return a[0] * (exp(-a[1] * x) + exp(-a[1] * (10.0 - x))) * sin(a[2] * x);
}

int main()
{
    float x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    mdp_measure y[10];
    int i, ma = 3;
    float a[3] = {1000.0, 0.1, 0.1};

    mdp_matrix covar(ma, ma);
    for (i = 0; i < 10; i++)
    {
        x[i] = i;
        y[i].set(f(x[i], a, ma, 0), 0.1);
    }

    a[0] = 1;
    a[1] = 1.0, a[2] = 0.2;

    BaesyanLevenbergMarquardt(x, y, 0, 10, a, ma, covar, f);

    for (i = 0; i < 10; i++)
        printf("%f %f (%f) %f\n",
               x[i], y[i].mean, y[i].error, f(x[i], a, ma, 0));

    printf("%f %f %f\n", a[0], a[1], a[2]);
    std::cout << covar;
}