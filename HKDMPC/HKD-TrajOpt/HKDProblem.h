#ifndef HKDPROBLEM_H
#define HKDPROBLEM_H

#include "TrajectoryManagement.h"
#include "SinglePhase.h"
#include "ConstraintsBase.h"
#include "HKDModel.h"
#include "HKDReset.h"
#include <memory>
#include <deque>
#include "QuadReference.h"
#include "HKDReference.h"

#include <lcm/lcm-cpp.hpp>
#include "hkd_problem_data_lcm_t.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

struct HKDPlanConfig
{
    float plan_duration;            // planning horizon in seconds
    float timeStep;                 // simulation timestep
    int nsteps_between_mpc;         // number of simulation time steps between mpc update (number of mpc steps control is applied)
};

template<typename T>
struct HKDProblemData
{
    QuadReference* quad_ref_ptr = nullptr;

    std::deque<shared_ptr<Trajectory<T,24,24,0>>> trajectory_ptrs;
    std::deque<shared_ptr<SinglePhase<T,24,24,0>>> phase_ptrs;   

    std::deque<int> phase_horizons;
    std::deque<bool> is_phase_reach_end;
    std::deque<float> phase_start_times;
    std::deque<float> phase_end_times;
    std::deque<VecM<double, 4>> contact_durations;
    std::deque<VecM<int, 4>> phase_contacts;
    int n_phases;

    void clear(){
        trajectory_ptrs.clear();
        phase_ptrs.clear();
        phase_horizons.clear();
        is_phase_reach_end.clear();
        phase_start_times.clear();
        phase_end_times.clear();
        phase_contacts.clear();
        contact_durations.clear();
        n_phases = 0;
        quad_ref_ptr = nullptr;
    }

    void pop_front_phase(){
        trajectory_ptrs.pop_front();
        phase_ptrs.pop_front();
        phase_start_times.pop_front();
        phase_end_times.pop_front();
        phase_horizons.pop_front();
        is_phase_reach_end.pop_front();
        phase_contacts.pop_front();
        contact_durations.pop_front();
        n_phases --;
    }
};

template<typename T>
inline void loadConstrintParameters(const std::string&fileName, 
                                    REB_Param_Struct<T>& GRF_reb_param,
                                    REB_Param_Struct<T>& Swing_reb_param,
                                    AL_Param_Struct<T>& TD_al_param)
{
	boost::property_tree::ptree pt;
    boost::property_tree::read_info(fileName, pt);
	std::cout << "********* loading MHPC Constraint Parameter from file **********\n" << fileName << "\n\n";

	GRF_reb_param.delta = pt.get<T>("GRF_ReB.delta");
	GRF_reb_param.delta_min = pt.get<T>("GRF_ReB.delta_min");
	GRF_reb_param.eps = pt.get<T>("GRF_ReB.eps");

    Swing_reb_param.delta = pt.get<T>("Swing_ReB.delta");
	Swing_reb_param.delta_min = pt.get<T>("Swing_ReB.delta_min");
	Swing_reb_param.eps = pt.get<T>("Swing_ReB.eps");	

    TD_al_param.sigma = pt.get<T>("TD_AL.sigma");
	TD_al_param.lambda = pt.get<T>("TD_AL.lambda");	
    TD_al_param.sigma_max = pt.get<T>("TD_AL.sigma_max");
}

template<typename T>
class HKDProblem
{
public:
    const static size_t xs = 24;
    const static size_t us = 24;
    const static size_t ys = 0;

public:
    HKDProblem(){plan_duration = 0; dt_sim = 0; nsteps_between_mpc = 0; dt_mpc = 0;}   

    void set_problem_data(HKDProblemData<T>* pdata_in, const HKDPlanConfig& config){
        plan_duration = config.plan_duration;
        dt_sim = config.timeStep;
        nsteps_between_mpc = config.nsteps_between_mpc;
        dt_mpc = dt_sim * nsteps_between_mpc;    

        pdata = pdata_in;
        quad_ref_ptr = pdata_in->quad_ref_ptr;                
    }

    void initialization();

    void update();    

    void clear_problem_data(){
        if (pdata != nullptr)
        {
            pdata->clear();
        }
            
    }    

    void pretty_print();

    void lcm_publish();

    void reset_lcm_data();

private:

    void create_problem_one_phase(shared_ptr<SinglePhase<T,24,24,0>> phase, int phase_idx);

    void add_tconstr_one_phase(shared_ptr<SinglePhase<T,24,24,0>> phase, int phase_idx);

private:
    HKDProblemData<T>* pdata = nullptr;    
    HKDSinglePhaseReference hkd_reference;
    HKD::Model<T> hkdModel;
    HKDReset<T> hkdReset;

    QuadReference* quad_ref_ptr = nullptr;

    // Define parameters for constraint management
    REB_Param_Struct<T> grf_reb_param;
    REB_Param_Struct<T> swing_reb_param;
    AL_Param_Struct<T> td_al_param;

    float plan_duration; 
    float dt_sim;                   // integratiion timestep
    int nsteps_between_mpc;     // number of integration timesteps between mpc updates    
    float dt_mpc;                   // time period between two mpc updates
    

    lcm::LCM _lcm;
    hkd_problem_data_lcm_t lcm_pdata;
};

#endif