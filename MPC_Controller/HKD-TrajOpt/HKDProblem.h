#ifndef HKDPROBLEM_H
#define HKDPROBLEM_H

#include "HKDContactSchedule.h"
#include "TrajectoryManagement.h"
#include "SinglePhase.h"
#include "ConstraintsBase.h"
#include "CostBase.h"
#include "HKDModel.h"
#include "HKDReset.h"
#include <memory>
#include <deque>
#include "HKDReference.h"

#include <lcm/lcm-cpp.hpp>
#include "hkd_problem_data_lcm_t.hpp"

using std::deque;

template<typename T>
struct HKDPlanConfig
{
    T plan_duration; 
    T timeStep;                 // integratiion timestep
    int nsteps_between_mpc;     // number of time steps between mpc update (number of mpc steps control is applied)
};

template<typename T>
struct HKDProblemData
{
    HKDReference<T>* reference_ptr;
    HKDReferenceData<T>* ref_data_ptr;

    deque<shared_ptr<Trajectory<T,24,24,0>>> trajectory_ptrs;
    deque<shared_ptr<SinglePhase<T,24,24,0>>> phase_ptrs;   
};


template<typename T>
class HKDProblem
{
public:
    const static size_t xs = 24;
    const static size_t us = 24;
    const static size_t ys = 0;

private:
    HKDProblemData<T>* pdata;
    HKDReference<T>* reference;
    HKDReferenceData<T>* ref_data;

    HKDModel<T> hkdModel;
    HKDReset<T> hkdReset;

    // initialize parameters for constraint management
    REB_Param_Struct<T> grf_reb_param;
    REB_Param_Struct<T> swing_reb_param;
    AL_Param_Struct<T> td_al_param;

    deque<deque<VecM<T,xs>>> Xbar_prev;
    deque<deque<VecM<T,us>>> Ubar_prev;

    T plan_duration; 
    T timeStep;             // integratiion timestep
    T time_between_mpc;     // time period between two mpc updates
    int nsteps_between_mpc;  // number of time steps between mpc updates

    lcm::LCM _lcm;
    hkd_problem_data_lcm_t lcm_pdata;

public:
    HKDProblem(){}
    HKDProblem(HKDProblemData<T>* pdata_in, const HKDPlanConfig<T>& config) {
        setup(pdata_in, config);
    }

    void setup(HKDProblemData<T>* pdata_in, const HKDPlanConfig<T>& config){
        plan_duration = config.plan_duration;
        timeStep = config.timeStep;
        nsteps_between_mpc = config.nsteps_between_mpc;
        time_between_mpc = timeStep * nsteps_between_mpc;    

        pdata = pdata_in;
        reference = pdata->reference_ptr;
        ref_data = pdata->ref_data_ptr;
    }

    void initialization();

    void update();    

    void create_problem_one_phase(shared_ptr<SinglePhase<T,24,24,0>> phase, int idx, bool add_tconstr = true);

    void add_tconstr_one_phase(shared_ptr<SinglePhase<T,24,24,0>> phase, int idx);

    void print();

    void lcm_publish();

    void reset_lcm_data();
};

#endif