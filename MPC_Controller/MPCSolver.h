#ifndef OPT_MAIN_H
#define OPT_MAIN_H

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
#include "Imitation_Reference.h"
#include "HKDReference.h"

template<typename T>
class MPCSolver
{
public:
    MPCSolver() : mpc_lcm(getLcmUrl(255))
    {
        // Setup reference
        string imitation_path = "../PolicyRollout/";
        string contact_fname = imitation_path + "contact_post.csv";
        string state_fname = imitation_path + "state_post.csv";
        imitation_ref.load_contact_data(contact_fname);
        imitation_ref.load_state_data(state_fname);
        imitation_ref.compute_status_duration();

        opt_ref.set_topLevel_reference(imitation_ref.get_data_ptr());

        opt_problem_data.reference_ptr = &opt_ref;
        opt_problem_data.ref_data_ptr = opt_ref.get_referenceData_ptr();

        // Default DDP options
        ddp_options.update_penalty = 5;
        ddp_options.update_relax = 1;
        ddp_options.update_ReB = 1;
        ddp_options.update_regularization = 4;      
        ddp_options.cost_thresh = 1e-02;
        ddp_options.AL_active = 1;
        ddp_options.ReB_active = 0;
        ddp_options.pconstr_thresh = .003;
        ddp_options.tconstr_thresh = .003;
        ddp_options.MS = true;
        ddp_options.merit_rho = 5*1e02;


        // Check LCM initialization
        if (!mpc_lcm.good())
        {
            printf(RED);
            printf("Failed to inialize mpc_lcm for hkd command\n");
            return;
        }       

        mpc_lcm.subscribe("mpc_data", &MPCSolver::mpcdata_lcm_handler, this);       
        solve_time = 0;
    }
    void mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                             const hkd_data_lcmt *msg);
    void publish_mpc_cmd();
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
    HKDReference<T> opt_ref;
    Imitation_Reference<T> imitation_ref; 

    HKDPlanConfig<T> mpc_config;
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

    // foot placement
    Vec3<float> pf[4];

    // mutex lock
    std::mutex mpc_mutex;    

    // solve time
    float solve_time;
};



#endif