## **Dependency**
This implementation uses [Eigen](https://gitlab.com/libeigen/eigen) for linear algegra, [Pinocchio](git@github.com:stack-of-tasks/pinocchio.git) for analytical derivatives, [LCM](https://github.com/lcm-proj/lcm/releases) for communications to low-level controllers, and [boost](https://www.boost.org/users/history/) for reading configuration files. A customized **Hybrid-Systems DDP (HS-DDP)** solver is employed to solve the nonlinear trajectory optimization problem. A C++ implementation of **HS-DDP** is included in this repo and can be found [here](https://github.com/heli-sudoo/HKD-MPC/tree/ICRA22%2BIROS23/MPC_Controller/HSDDPSolver).

- [Eigen3](https://gitlab.com/libeigen/eigen)
- [LCM1.4.0](https://github.com/lcm-proj/lcm/releases)
- [Boost1.71](https://www.boost.org/users/history/)
- [Pinocchio](git@github.com:stack-of-tasks/pinocchio.git)


## **Build**
Once Eigen and LCM are successfully installed, generate necessary lcm types

```bash
cd scripts
./make_types.sh
```

To build the MPC controller

```bash
mkdir build && cd build
cmake ..
make -j4
```

To run the HKDMPC controll
```bash
cd build
HKDMPC/mpc_solve
```