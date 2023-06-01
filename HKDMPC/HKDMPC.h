#ifndef HKDMPC_H
#define HKDMPC_H

#include <thread>
#include <mutex>
#include <lcm/lcm-cpp.hpp>
#include "hkd_command_lcmt.hpp"
#include "hkd_data_lcmt.hpp"
#include "opt_sol_lcmt.hpp"

#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"
#include "HKDModel.h"
#include "HKDProblem.h"
#include "HKDReset.h"
#include "MultiPhaseDDP.h"
#include "cTypes.h"
#include "utilities.h"
#include "HKDReference.h"
#include "solver_info_lcmt.hpp"

template<typename T>
class HKDMPCSolver
{
public:
    HKDMPCSolver() : mpc_lcm(getLcmUrl(255)),
                  solver_info_lcm(getLcmUrl(255))
    {
        // Setup reference
        string reference_file_path = "../Reference/Data/quad_reference.csv";       

        quad_reference.load_top_level_data(reference_file_path);                        

        // Check LCM initialization
        if (!mpc_lcm.good())
        {
            printf(RED);
            printf("Failed to inialize mpc_lcm for hkd command\n");
            return;
        }       

        mpc_lcm.subscribe("mpc_data", &HKDMPCSolver::mpcdata_lcm_handler, this);       
        solve_time = 0;
    }
    void mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                             const hkd_data_lcmt *msg);
    void publish_mpc_cmd();
    void publish_solver_info();
    void publish_debugfoot();
    void initialize();
    void update();
    void update_foot_placement();
    void run(){
        while (mpc_lcm.handle()==0){}
    }

public:
    // MPC
    HKDProblem<T> opt_problem;
    HKDProblemData<T> opt_problem_data;
    QuadReference quad_reference;

    HKDPlanConfig mpc_config;
    HSDDP_OPTION ddp_options;

    T dt_mpc;
    T mpc_time;
    T mpc_time_prev;
    int mpc_iter;
    
    DVec<T> xinit;
    VecM<T, 12> body, qdummy, qJ;
    VecM<T, 3> pos, eul, vel, omega;

    // LCM message
    hkd_data_lcmt hkd_data;
    hkd_command_lcmt hkd_cmds;
    opt_sol_lcmt debug_foot_data;
    lcm::LCM mpc_lcm;
    lcm::LCM solver_info_lcm;
    solver_info_lcmt solver_info;

    // foot placement
    Vec3<float> pf[4];

    // mutex lock
    std::mutex mpc_mutex;    

    // solve time
    float solve_time;
};



#endif