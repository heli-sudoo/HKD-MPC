#include "HKDMPC.h"

int main()
{
    HKDMPCSolver<double> mpc;    
    mpc.initialize();
    mpc.run();
}