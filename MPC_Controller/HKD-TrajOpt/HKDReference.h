#ifndef HKDREFERENCE_H
#define HKDREFERENCE_H

#include <deque>
#include <memory>
#include <algorithm>
#include <cassert>
#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"

using std::deque;
using std::shared_ptr;

template<typename T>
struct HKDReferenceData
{
    deque<VecM<int, 4>> contactSeq; // should have at least one more element than n_phases
    DMat<T> statusDuration;
    deque<T> startTimes;    // start times of each contact phase
    deque<T> endTimes;      // end times of each contact phase
    deque<int> horizons;    // horizon of each contact phase
    T dt;
    deque<deque<VecM<T, 24>>> Xr;
    deque<deque<VecM<T, 24>>> Ur;
    deque<deque<VecM<T, 0>>> Yr;
    int n_phases;

    void clear()
    {
        contactSeq.clear();        
        startTimes.clear();
        endTimes.clear();
        horizons.clear();
        Xr.clear();
        Ur.clear();
        Yr.clear();
        dt = 0;
        n_phases = 0;
    }

    void resize(int n_phases_in)
    {
        n_phases = n_phases_in;
        contactSeq.resize(n_phases+1);
        statusDuration.resize(4, n_phases+1);
        startTimes.resize(n_phases);
        endTimes.resize(n_phases);
        horizons.resize(n_phases);
        Xr.resize(n_phases);
        Ur.resize(n_phases);
        Yr.resize(n_phases);
    }
    int size()
    {
        return n_phases;
    }
    void pop_front()
    {
        contactSeq.pop_front();
        startTimes.pop_front();
        endTimes.pop_front();
        horizons.pop_front();
        Xr.pop_front();
        Ur.pop_front();
        Yr.pop_front();
        n_phases--;
    }
};


template<typename T>
class HKDReference
{
private:
    HKDReferenceData<T> data;                   // reference data used in optimization
    HKDReferenceData<T> *tp_data_ptr = nullptr;    // reference data of long horizon from top level planning (or offline)
    T dt;
    int last_phase_idx;     // keep track the index of the last phase w.r.t. the top level data

public:
    HKDReference(/* args */){data.clear();dt = 0;}

    T get_lastphase_endtime_tp() {return tp_data_ptr->endTimes[last_phase_idx];}
    
    void set_topLevel_reference(HKDReferenceData<T> *tp_data_ptr_in){
        tp_data_ptr = tp_data_ptr_in;
        dt = tp_data_ptr->dt;
        }

    HKDReferenceData<T> * get_referenceData_ptr() {return &data;}

    void clear() {data.clear(); last_phase_idx = 0;}

    void initialize_referenceData(T duration){
        // induce the number of phases and the horizon of each phase
        data.clear();
        data.dt = dt;
        deque<int>& horizons = data.horizons;
        int i(0);
        long int duration_ns = 1e9*duration; // conver to nano second
        long int timeStep_ns = 1e9*dt;
        while (i < tp_data_ptr->n_phases)
        {
            T startTime_i = tp_data_ptr->startTimes[i];
            T endTime_i = tp_data_ptr->endTimes[i];
            long int startTime_i_ns = 1e9*startTime_i;
            long int endTime_i_ns = 1e9*endTime_i;
            if (duration_ns == endTime_i_ns)
            {
                data.startTimes.push_back(startTime_i);
                data.endTimes.push_back(endTime_i);
                horizons.push_back(tp_data_ptr->horizons[i]);
                break;
            }
            
            if (duration_ns < endTime_i_ns)
            {
                data.startTimes.push_back(startTime_i);
                data.endTimes.push_back(duration);
                horizons.push_back((duration_ns-startTime_i_ns)/timeStep_ns);
                break;               
            }
            data.startTimes.push_back(startTime_i);
            data.endTimes.push_back(endTime_i);
            horizons.push_back(tp_data_ptr->horizons[i]);
            i++;            
        }
        
        int n_phases = horizons.size();
        int iter_offset = n_phases+1;
        data.resize(n_phases);
        std::copy(tp_data_ptr->contactSeq.begin(), tp_data_ptr->contactSeq.begin()+iter_offset,
                data.contactSeq.begin());
        data.statusDuration = tp_data_ptr->statusDuration.leftCols(iter_offset);

        // For each phase, copy the data of proper length from top-level reference info to
        // the local reference info being used in optimization
        for (int i = 0; i < n_phases; i++)
        {
            data.endTimes[i] = data.startTimes[i] + dt * horizons[i];
            std::copy(tp_data_ptr->Xr[i].begin(), tp_data_ptr->Xr[i].begin()+horizons[i]+1,
                    std::back_inserter(data.Xr[i]));
            std::copy(tp_data_ptr->Ur[i].begin(), tp_data_ptr->Ur[i].begin()+horizons[i],
                    std::back_inserter(data.Ur[i]));
            std::copy(tp_data_ptr->Yr[i].begin(), tp_data_ptr->Yr[i].begin()+horizons[i],
                    std::back_inserter(data.Yr[i]));
        }
        data.n_phases = n_phases;
        last_phase_idx = n_phases - 1;

        // compute_status_duration();
    }

    // Shift the reference data forward by one time step
    void step(){
        // If the first phase reaches to the end, pop front the first phase
        if (data.horizons.front()<=1){
            data.pop_front();
        }else{ //else keep the first phase but reduce the state, control, and output
            // printf("Xr[0] size = %lu \n", data.Xr[0].size());
            // printf("starttime = %f \n", data.startTimes.front());
            // printf("horizon = %i \n", data.horizons.front());
            data.Xr[0].pop_front();
            data.Ur[0].pop_front();
            data.Yr[0].pop_front();
            data.startTimes.front()+=dt;
            data.horizons.front()--;
        }
        
        // If the last phase goes beyond the end time, add one more phase
        auto tp_ptr = tp_data_ptr;
        if (data.endTimes.back()+dt > tp_ptr->endTimes[last_phase_idx]+1e-08)
        {
            data.resize(data.size()+1);
            data.contactSeq.back()=tp_ptr->contactSeq[last_phase_idx+2];
            data.startTimes.back()=tp_ptr->startTimes[last_phase_idx+1];
            data.endTimes.back() = data.startTimes.back() + dt;
            data.horizons.back() = 1;

            std::copy(tp_ptr->Xr[last_phase_idx+1].begin(), tp_ptr->Xr[last_phase_idx+1].begin()+2,
                std::back_inserter(data.Xr.back()));
            data.Ur.back().push_back(tp_ptr->Ur[last_phase_idx+1].front());
            data.Yr.back().push_back(tp_ptr->Yr[last_phase_idx+1].front());

            last_phase_idx++;
        }else{                        
            data.Xr.back().push_back(tp_ptr->Xr[last_phase_idx][data.horizons.back()+1]);
            data.Ur.back().push_back(tp_ptr->Ur[last_phase_idx][data.horizons.back()]);
            data.Yr.back().push_back(tp_ptr->Yr[last_phase_idx][data.horizons.back()]);

            data.endTimes.back()+=dt;
            data.horizons.back()++;
        }
        data.statusDuration = tp_data_ptr->statusDuration.middleCols(last_phase_idx-(data.n_phases-1), data.n_phases);        
    }
};




#endif //HKDREFERENCE_H