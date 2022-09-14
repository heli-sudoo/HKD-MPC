#include "MPCSolver.h"

int main()
{
    MPCSolver<double> mpc;
    mpc.initialize();
    mpc.run();
}