## **HKD-MPC**
The HKD-MPC is a nonlinear MPC controller for agile and versatile quadruped locomotion. It enables highly dynamic jumping and rapid gait transition in one framework, without the need to switch between multiple controllers. The HKD-MPC has been tested on Unitree A1 and MIT Mini Cheetah. This repo contains the 
implementation for the MIT Mini Cheetah. 

## **Demo**
<img src="demo/A1_jump.gif" height="130">  <img src="demo/A1_hop.gif" height="130"> <img src="demo/MC_hop.gif" height="130">

## **Dependency**
This implementation only requires two external dependencies, Eigen and LCM. A customized **Hybrid-Systems DDP (HS-DDP)** solver is employed to solve the nonlinear trajectory optimization problem. A C++ implementation of **HS-DDP** is included in this repo and can be found [here](https://github.com/heli-sudoo/HKD-MPC/tree/ICRA22%2BIROS23/MPC_Controller/HSDDPSolver).
- [Eigen3](https://gitlab.com/libeigen/eigen)
- [LCM1.4.0](https://github.com/lcm-proj/lcm/releases)


## **Build**
Once Eigen and LCM are successfully installed, generate necessary lcm types

```bash
cd scripts
./make_types.sh
```

Then build the controller

```bash
cd ../MPC_Controller
mkdir build && cd build
cmake ..
make -j4
```

## **Run in Cheetah Software**
Instructions coming soon

## **Publications**
[1] Li, H., Yu, W., Zhang, T., & Wensing, P. M. (2022, October). Zero-shot retargeting of learned quadruped locomotion policies using hybrid kinodynamic model predictive control. In 2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS) (pp. 11971-11977). IEEE.

[2] Li, H., Zhang, T., Yu, W., & Wensing, P. M. (2022). Versatile Real-Time Motion Synthesis via Kino-Dynamic MPC with Hybrid-Systems DDP. arXiv preprint arXiv:2209.14138.