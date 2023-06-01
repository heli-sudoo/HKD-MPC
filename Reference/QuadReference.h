#pragma once
#ifndef QUADRUPED_REFERENCE_H
#define QUADRUPED_REFERENCE_H

#include <deque>
#include <memory>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"

struct QuadAugmentedState
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VecM<double, 12> body_state;      // [eul, pos, ang, vel]
    VecM<double, 12> qJ;              // joint angle
    VecM<double, 12> qJd;             // joint velocity
    VecM<double, 12> foot_placements; // foothold position reference (zero if in swing) in world frame
    VecM<double, 12> grf;             // GRF reference (zero if in swing) in world frame
    VecM<double, 12> torque;          // Joint torque reference
    VecM<int, 4> contact;            // contact status
    VecM<double, 4> status_dur;       // contact status duration

    void SetZero()
    {
        body_state.setZero();
        qJ.setZero();
        qJd.setZero();
        foot_placements.setZero();
        grf.setZero();
        torque.setZero();
        contact.setZero();
        status_dur.setZero();
    }

    void print()
    {
        std::cout << "body state = \n";
        std::cout << body_state.transpose() << "\n";

        std::cout << "qJ = \n";
        std::cout << qJ.transpose() << "\n";

        std::cout << "qJd = \n";
        std::cout << qJd.transpose() << "\n";

        std::cout << "foot_placements = \n";
        std::cout << foot_placements.transpose() << "\n";

        std::cout << "grf = \n";
        std::cout << grf.transpose() << "\n";

        std::cout << "torque = \n";
        std::cout << torque.transpose() << "\n";

        std::cout << "contact = \n";
        std::cout << contact.transpose() << "\n";

        std::cout << "status_dur = \n";
        std::cout << status_dur.transpose() << "\n";
    }
};

class QuadReferenceData
{
public:
    QuadReferenceData()
    {
        dt = 0;
        start_time = 0;
        end_time = 0;
    } // Default constructor: create empty data

    QuadReferenceData(size_t sz) : astates(sz)
    {
        dt = 0;
        start_time = 0;
        end_time = 0;
    }

    float get_duration() { return end_time - start_time; } // get the duration of the trajectory

    size_t size() const { return astates.size(); }                  // get the size of the trajectory

    void clear()
    {
        astates.clear();
    }

    void pop_front()
    {
        astates.pop_front();
    }

    void push_back(QuadAugmentedState& astate_to_add)
    {
        astates.push_back(astate_to_add);
    }

    QuadAugmentedState& operator[](size_t k)
    {
        assert((k < astates.size()));
        return astates[k];
    }    

    QuadAugmentedState& at(size_t k)
    {
        assert((k < astates.size()));
        return astates[k];
    }

    QuadAugmentedState* get_ptr(size_t k)
    {
        return &astates[k];
    }

    std::deque<QuadAugmentedState>::iterator begin()
    {
        return astates.begin();
    }

    std::deque<QuadAugmentedState>::iterator end()
    {
        return astates.end();
    }

    std::deque<QuadAugmentedState>& get_container()
    {
        return astates;
    }

public:
    float dt; // discretization time step of the reference data
    float start_time;
    float end_time;

private:
    std::deque<QuadAugmentedState> astates;
};

class QuadReference
{
private:
    QuadReferenceData data;    // reference data used in optimization
    QuadReferenceData tp_data; // very-long horizon reference data (currently loaded from external files)
    float t_cur;               // track where we are in the top-level reference (second)
    int k_cur;                 // track where we are in the top level reference (index)
    float dur;                 // overall planning horizon
    float dt;                  // timestep used to discretize the reference
    int sz;                    // size of the reference trajectory (dur/dt + 1) (i.e., number of points)

public:
    QuadReference()
    {
        dt = 0;
        dur = 0;
        sz = 0;
    }

    void load_top_level_data(const std::string &fname, bool reorder = false);

    QuadAugmentedState* get_a_reference_ptr_at_t(float t);

    QuadReferenceData* const get_data_ptr() {return &data;}

    QuadReferenceData* const get_tp_data_ptr() {return &tp_data;}

    void initialize(float plan_horizon);

    void step(float dt_sim);

    void update_foot_placements(const VecM<double, 12>& body_state_t);

    void get_contact_at_t(VecM<int, 4> &contact, float t);

    void get_contact_duration_at_t(VecM<double, 4> &contact_dur, float t);
   
    int get_data_size() { return sz; } // Get the length of the reference trajectory

    float get_dt() { return dt; } // Get dt in reference trajectory

    float get_start_time() { return data.start_time; } // Start time of the reference relative to the long-horizon top-level reference data

    float get_end_time() { return data.end_time; } // End time of the reference relative to the long-horizon top-level reference data

private:
    void reorder_states();
};


#endif // QUADRUPED_REFERENCE_H