#pragma once
#ifndef HSDDP_COMPOUNDTYPES_H
#define HSDDP_COMPOUNDTYPES_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include "HSDDP_CPPTypes.h"
#include <iostream>
#include <cstring>


#ifdef TIME_BENCHMARK
extern std::vector<TIME_PER_ITERATION> time_ddp;
#endif //TIME_BENCHMARK


//  HSDDP parameters
struct HSDDP_OPTION
{
    double alpha = 0.1;               // line search udpate paramer
    double gamma = 0.1;              // scale the expected cost reduction
    double update_penalty = 8;        // penalty update parameter
    double update_relax = 0.1;        // relaxation parameter udpate parameter
    double update_regularization = 2; // regularization parameter update parameter
    double update_ReB = 7;            // update barrier function weighting
    int max_DDP_iter = 3;             // maximum inner loop iteration (for the initial solve)
    int max_AL_iter = 2;              // maximum outer loop iteration (for the initial solve)
    int max_DDP_iter_runtime = 1;     // maximum inner loop iteration (in MPC loop)
    int max_AL_iter_runtime = 2;      // maximum outer loop iteration (in MPC loop)
    double cost_thresh = 1e-03;        // inner loop convergence threshhold
    double tconstr_thresh = 1e-03;
    double pconstr_thresh = 1e-03;
    double dynamics_feas_thresh = 1e-03;    // threshold to accept dynamics infeasibility
    double merit_rho = 1e04;          // weighting parameter for merit function
    double merit_scale = 0.2;
    double merit_offset = 10;
    bool AL_active = 1;               // activate terminal constraint
    bool ReB_active = 1;              // activate path constraint
    bool smooth_active = 0;           // activate control smoothness penalization
    bool MS = true;                   // shooting method (true MS, false SS)
    int nsteps_per_node = 1;          // define a shooting node for every nsteps_per_node integration time steps

    void print(){
        std::cout << "===================== HSDDP Setting =================== \n";
        std::cout << "Multiple Shooting \t" << MS << "\n";
        std::cout << "Number of integration time steps per node \t" << nsteps_per_node << "\n";
        std::cout << "Is Terminal Constraint active \t" << AL_active << "\n";
        std::cout << "Is path constraint active \t" << ReB_active << "\n";
        std::cout << "Maximum inner-loop iterations \t" << max_DDP_iter << "\n";
        std::cout << "Maximum outer-loop iterations \t" << max_AL_iter << "\n";
        std::cout << "Line search update param \t" << alpha << "\n";
        std::cout << "Merit function penalty \t" << merit_rho << "\n";
        std::cout << "Merit scale parameter \t" << merit_scale << "\n";
        std::cout << "Merit offset parameter \t" << merit_offset << "\n";
        std::cout << "Terminal constraint threshold \t" << tconstr_thresh << "\n";
        std::cout << "Path constraint threshold \t" << pconstr_thresh << "\n";
        std::cout << "Dynamics infeasibility threshold \t" << dynamics_feas_thresh << "\n";
        std::cout << "Cost convergence threshold \t" << cost_thresh << "\n\n";
    }
};

inline void loadHSDDPSetting(const std::string& filename, HSDDP_OPTION& setting){
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);

    std::cout << "********* loading HSDDP setting from file **********\n" << filename << "\n\n";
    setting.alpha = pt.get<double>("ddp.alpha");
    setting.gamma = pt.get<double>("ddp.gamma");
    setting.update_penalty = pt.get<double>("ddp.update_penalty");
    setting.update_relax = pt.get<double>("ddp.update_relax");
    setting.update_ReB = pt.get<double>("ddp.update_ReB");
    setting.max_DDP_iter = pt.get<int>("ddp.max_DDP_iter");
    setting.max_AL_iter = pt.get<int>("ddp.max_AL_iter");
    setting.max_DDP_iter_runtime = pt.get<int>("ddp.max_DDP_iter_runtime");
    setting.max_AL_iter_runtime = pt.get<int>("ddp.max_AL_iter_runtime");
    setting.cost_thresh = pt.get<double>("ddp.cost_thresh");
    setting.tconstr_thresh = pt.get<double>("ddp.tconstr_thresh");
    setting.pconstr_thresh = pt.get<double>("ddp.pconstr_thresh");
    setting.dynamics_feas_thresh = pt.get<double>("ddp.dynamics_feas_thresh");
    setting.merit_rho = pt.get<double>("ddp.merit_rho");
    setting.merit_scale = pt.get<double>("ddp.merit_scale");
    setting.merit_offset = pt.get<double>("ddp.merit_offset");
    setting.AL_active = pt.get<bool>("ddp.AL_active");
    setting.ReB_active = pt.get<bool>("ddp.ReB_active");
    setting.MS = pt.get<bool>("ddp.MS");
    setting.nsteps_per_node = pt.get<int>("ddp.nsteps_per_node");
}

/* Running Cost Data Structure */
template<typename T, size_t xs_, size_t us_, size_t ys_>
struct RCostData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T l;
    VecM<T, xs_> lx;
    VecM<T, us_> lu;
    VecM<T, ys_> ly;
    MatMN<T, xs_, xs_> lxx;
    MatMN<T, us_, xs_> lux;
    MatMN<T, us_, us_> luu;
    MatMN<T, ys_, ys_> lyy;

    RCostData(){Zeros();}
    void Zeros()
    {
        l = 0;
        lx.setZero();
        lu.setZero();
        ly.setZero();
        lxx.setZero();
        luu.setZero();
        lux.setZero();
        lyy.setZero();
    }
    void add(const RCostData<T,xs_,us_,ys_>& cost_to_add)
    {
        l += cost_to_add.l;
        lx += cost_to_add.lx;
        lu += cost_to_add.lu;
        ly += cost_to_add.ly;
        lxx += cost_to_add.lxx;
        luu += cost_to_add.luu;
        lux += cost_to_add.lux;
        lyy += cost_to_add.lyy;
    }
};

/* Terminal Cost Data Structure */
template<typename T, size_t xs_>
struct TCostData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T Phi;
    VecM<T, xs_> Phix;
    MatMN<T, xs_, xs_> Phixx;

    TCostData(){Zeros();}
    void Zeros()
    {
        Phi = 0;
        Phix.setZero();
        Phixx.setZero();
    }
    void add(const TCostData<T,xs_>& cost_to_add)
    {
        Phi += cost_to_add.Phi;
        Phix += cost_to_add.Phix;
        Phixx += cost_to_add.Phixx;
    }
};


#endif // HSDDP_COMPOUNDTYPES_H