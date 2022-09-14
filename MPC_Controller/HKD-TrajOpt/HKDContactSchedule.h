#ifndef HKD_CONTACT_SCHEDULE_H
#define HKD_CONTACT_SCHEDULE_H

#include <vector>
#include <deque>
#include "HSDDP_CPPTypes.h"
#include "TrajectoryManagement.h"

using std::vector;
using std::deque;

template<class Model_>
class ContactSchedule
{
private:
    typedef typename Model_::Scalar valType;
    static const int xs = Model_::state_dim;
    static const int us = Model_::contrl_dim;
    static const int ys = Model_::output_dim;

public:
    ContactSchedule(){}
    void set_timestep(valType timestep_) {dt = timestep_;}
    void get_phase_indices(vector<size_t>& indices, valType t, valType plan_duration);
    void get_phase_durations(deque<valType>& phase_durations, valType t, valType plan_duration);
    void get_phase_horizons(deque<size_t>&phase_horizons, valType t, valType plan_duration, valType timestep);
    void get_phase_contacts(vector<VecM<size_t,4>>&phase_ctacts, valType t, valType plan_duration);
    void get_phase_reference_it(vector<BriefTrajectoryIterator<Model_>>&ref_iters, valType t, valType plan_duration);
    void get_phase_horizons(deque<size_t>&horizons_new, valType told, valType tnew, valType plan_dur_old, valType plan_dur_new, valType timestep);
    void load_contact_data(const string&); 
    void load_state_data(const string&);
    void update() {} // Placeholder
    void print_info();
    void print_state();
    void log_reference();
    size_t get_contactStatus(valType t_cur, int foot);
    valType get_statusDuration(valType t_cur, int foot);
    
public:
    vector<VecM<size_t, 4>> ctact_seq; // contact sequence
    vector<valType> startTime;
    vector<valType> endTime;
    vector<size_t> horizons;
    vector<BriefTrajectory<Model_>> contactTraj;    
    valType dt = 0.011;
    size_t num_phases_tot = 0;

    VecM<valType, 4> stanceTimes;
    VecM<valType, 4> swingTimes;

private:
    vector<size_t> startidx;
    vector<size_t> endidx;
};
#endif