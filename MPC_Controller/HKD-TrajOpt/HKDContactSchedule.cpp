#include "HKDContactSchedule.h"
#include <fstream>
#include <iostream>
#include "HSDDP_Utils.h"
#include "cTypes.h"
using Eigen::IOFormat;
using std::string;

template <class Model_>
void ContactSchedule<Model_>::load_contact_data(const string &ctact_fname)
{
    fstream fstrm(ctact_fname);
    string line; // get each row of the file to line
    if (!fstrm.is_open())
    {
        printf("Failed to open contact file \n");
        return;
    }
    // Get the first line and print the names of each column
    getline(fstrm, line);
    // Go through the rest of lines
    startTime.clear();
    endTime.clear();
    startidx.clear();
    endidx.clear();
    ctact_seq.clear();
    VecM<size_t, 4> contact;
    size_t sidx;
    size_t eidx;
    valType sT, eT;
    num_phases_tot = 0;
    string word;
    while (getline(fstrm, line))
    {
        stringstream s(line); // break line to s with delimeter (,)
        // std::cout << line << std::endl;
        int i = 0;
        while (getline(s, word, ','))
        {
            if (i <= 3)
            {
                contact[i] = stoul(word);
            }
            if (i == 4)
            {
                sT = static_cast<valType>(stod(word));
            }
            if (i == 5)
            {
                eT = static_cast<valType>(stod(word));
            }
            if (i == 6)
            {
                sidx = stoul(word);
            }
            if (i == 7)
            {
                eidx = stoul(word);
            }
            i++;
        }
        ctact_seq.push_back(contact);
        startidx.push_back(sidx);
        endidx.push_back(eidx);
        startTime.push_back(sT);
        endTime.push_back(eT);
        num_phases_tot++;
    }
    fstrm.close();
    horizons.clear();
    /* Process phase durations */
    for (size_t i = 0; i < ctact_seq.size(); i++)
    {
        horizons.push_back(endidx[i] - startidx[i] + 1);
    }
}

template <class Model_>
void ContactSchedule<Model_>::load_state_data(const string &state_fname)
{
    fstream fstrm(state_fname);
    if (!fstrm.is_open())
    {
        printf("Failed to open state file \n");
        return;
    }
    contactTraj.clear();
    string line;
    int line_count = 0;
    VecM<valType, xs> x;
    VecM<valType, us> u;
    VecM<valType, ys> y;
    y.setZero();
    VecM<valType, 3> grf;
    grf << 0.0, 0.0, 20.0;
    u << grf, grf, grf, grf, VecM<valType, 12>::Zero();
    BriefTrajectory<Model_> traj_one_phase;
    string word;
    for (size_t pidx = 0; pidx < startidx.size(); pidx++)
    {
        traj_one_phase.clear();
        while (line_count <= endidx[pidx])
        {
            getline(fstrm, line);
            stringstream s(line);
            size_t i = 0;
            while (getline(s, word, ','))
            {
                x[i] = static_cast<valType>(stod(word));
                i++;
            }
            traj_one_phase.X.push_back(x);
            traj_one_phase.Y.push_back(y);
            if (line_count < endidx[pidx])
            {
                traj_one_phase.U.push_back(u);
            }

            line_count++;
        }
        contactTraj.push_back(traj_one_phase);
    }
    fstrm.close();
}

/*
    @brief  Compute contact status for plan phases. The size of phase_ctacts = num_plan_phases + 1 to induce touch down status
    @params
            t: current time
            plan_duration: plan duration in seconds
*/
template <class Model_>
void ContactSchedule<Model_>::get_phase_contacts(vector<VecM<size_t, 4>> &phase_ctacts, valType t, valType plan_duration)
{
    phase_ctacts.clear();
    vector<size_t> phase_indices;
    get_phase_indices(phase_indices, t, plan_duration);
    for (auto &pidx : phase_indices)
    {
        phase_ctacts.push_back(ctact_seq[pidx]);
    }
    if (!phase_ctacts.empty())
    {
        phase_ctacts.push_back(ctact_seq[phase_indices.back() + 1]);
    }
}

template <class Model_>
void ContactSchedule<Model_>::get_phase_durations(deque<valType> &phase_durations, valType t, valType plan_duration)
{
    phase_durations.clear();
    vector<size_t> phase_indices;
    get_phase_indices(phase_indices, t, plan_duration);
    valType phase_dur;
    for (size_t pidx : phase_indices)
    {
        if (pidx == phase_indices.front())
        {
            phase_dur = endTime[pidx] - t;
        }
        else if (pidx == phase_indices.back())
        {
            phase_dur = t + plan_duration - startTime[pidx];
        }
        else
        {
            phase_dur = endTime[pidx] - startTime[pidx];
        }
        if (phase_dur+1e-8 >= dt)
        {
            phase_durations.push_back(phase_dur);
        }                
    }
}
template <class Model_>
void ContactSchedule<Model_>::get_phase_horizons(deque<size_t> &phase_horizons, valType t, valType plan_duration, valType timestep)
{
    phase_horizons.clear();
    deque<valType> phase_durations;
    get_phase_durations(phase_durations, t, plan_duration);
    for (size_t i = 0; i < phase_durations.size(); i++)
    {
        phase_horizons.push_back(size_t(round(phase_durations[i] / timestep)));
    }
}
template <class Model_>
void ContactSchedule<Model_>::get_phase_horizons(deque<size_t> &phase_horizons_new, valType told, valType tnew,
                                                 valType plan_dur_old, valType plan_dur_new, valType timestep)
{
    deque<size_t> phase_horizons_old;
    vector<size_t> phase_ids_old;
    vector<size_t> phase_ids_new;
    get_phase_horizons(phase_horizons_old, told, plan_dur_old, timestep);
    get_phase_horizons(phase_horizons_new, tnew, plan_dur_new, timestep);
    get_phase_indices(phase_ids_old, told, plan_dur_old);
    get_phase_indices(phase_ids_new, tnew, plan_dur_new);
    // Iterate through the indices of the old phases until find a match of the first phase index of the new phases
    for (size_t i = 0; i < phase_ids_old.size(); i++)
    {
        if (phase_ids_old[i] < phase_ids_new[0])
        {
            phase_horizons_new.push_front(0);
        }
        else
        {
            break;
        }
    }
}

/*
    @brief  Compute the indices of planning phase in the reference contact trajectory
    @params
            t: current time
            plan_duration: plan duration in seconds
    @return
            indices: index of the plan phases
*/
template <class Model_>
void ContactSchedule<Model_>::get_phase_indices(vector<size_t> &indices, valType t, valType plan_duration)
{
    indices.clear();
    for (size_t i = 0; i < startTime.size(); i++)
    {
        if (t + 1e-8 >= endTime[i]) // add a very small number to address the precision problem
            continue;
        if (t + plan_duration <= (startTime[i]+ 1e-8))
            break;
        indices.push_back(i);
    }
}

template <class Model_>
void ContactSchedule<Model_>::get_phase_reference_it(vector<BriefTrajectoryIterator<Model_>> &ref_iters, valType t, valType plan_duration)
{
    ref_iters.clear();
    // get phase indices
    vector<size_t> phase_indices;
    get_phase_indices(phase_indices, t, plan_duration);
    size_t pidx = phase_indices[0];
    // get the trajectory iterator for the first phase with horizon offset
    valType time_offset = t - startTime[pidx];
    size_t horizon_offset = round(time_offset / dt);
    BriefTrajectoryIterator<Model_> phase_ref_iter;
    contactTraj[pidx].get_iterator(phase_ref_iter, horizon_offset);

    ref_iters.push_back(phase_ref_iter);
    // get the trajectory iterators for the remaining trajectories (offset = 0)
    for (size_t i = 1; i < phase_indices.size(); i++)
    {
        pidx = phase_indices[i];
        contactTraj[pidx].get_iterator(phase_ref_iter, 0);
        ref_iters.push_back(phase_ref_iter);
    }
}
/*
    @brief  Get the contact status given the current time for one foot
*/
template <class Model_>
size_t ContactSchedule<Model_>::get_contactStatus(valType t_cur, int foot)
{
    if (t_cur > endTime.back())
    {
        printf("Reached the end of the reference motion \n");
        return 0;
    }
    
    for (size_t i = 0; i < startTime.size(); i++)
    {
        if (t_cur+1e-08>= startTime[i] && (t_cur+1e-08< endTime[i]))
        {
            return ctact_seq[i][foot];
        }
        // if (t_cur < endTime[i] && t_cur+1e-08>endTime[i])
        // {
        //     printf(RED);
        //     printf("t_cur = %.15f, endtime = %.15f\n", t_cur, endTime[i]);
        //     printf(RESET);
        // }
        
    }
}
/*
    @brief Get the duration of the contact status that is found using the current time
*/
template <class Model_>
typename Model_::Scalar ContactSchedule<Model_>::get_statusDuration(valType t_cur, int foot)
{
    size_t ctact_status;
    valType status_starttime;
    valType status_endtime;
    if (t_cur > endTime.back())
    {
        printf("Reached the end of the reference motion \n");
        return 0;
    }
    
    // iterate over all phases
    for (int i = 0; i < startTime.size(); i++)
    {
        // find the current contact status
        if (t_cur+1e-08>= startTime[i] && t_cur+1e-08< endTime[i])
        {
            ctact_status =  ctact_seq[i][foot];
            // find the start time for this status until status is changed
            for (int j = i;  j >= 0; j--)
            {
                if (ctact_status == ctact_seq[j][foot])
                {
                    status_starttime = startTime[j];
                }
                else{
                    break;
                }
            }
            // find the end time for this status until status is changed
            for (int j = i; j < endTime.size(); j++)
            {
                if (ctact_status == ctact_seq[j][foot])
                {
                    status_endtime = endTime[j];
                }
                else{
                    break;
                }
            }
        }
    }
    return status_endtime - status_starttime;
}

template <class Model_>
void ContactSchedule<Model_>::print_info()
{
    // Print contact information
    printf("\n");
    printf("Total number of contact phases = %lu\n", num_phases_tot);
    printf("%5s %4s %4s %4s %4s %10s %10s %10s %10s",
           "PIdx", "FR", "RL", "HR", "HL", "StartTime", "EndTime", "StartIndex", "EndIndex");
    printf("\n");
    for (size_t pidx = 0; pidx < num_phases_tot; pidx++)
    {
        printf("%5lu %4lu %4lu %4lu %4lu %10.3f %10.3f %10lu %10lu",
               pidx,
               ctact_seq[pidx][0], ctact_seq[pidx][1], ctact_seq[pidx][2], ctact_seq[pidx][3],
               startTime[pidx], endTime[pidx], startidx[pidx], endidx[pidx]);
        printf("\n");
    }
}

template <class Model_>
void ContactSchedule<Model_>::print_state()
{
    // Print eul and CoM pos information
    printf("\n");
    printf("%8s %8s %8s %8s %8s %8s",
           "yaw", "pitch", "row", "px", "py", "pz");
    printf("\n");
    for (size_t pidx = 0; pidx < num_phases_tot; pidx++)
    {
        auto &bt = contactTraj[pidx];
        for (size_t k = 0; k < bt.length(); k++)
        {
            auto &x = bt.X[k];
            printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f",
                   x[0], x[1], x[2], x[3], x[4], x[5]);
            printf("\n");
        }
    }
}

template <class Model_>
void ContactSchedule<Model_>::log_reference()
{
    string folder_name = "../Log/";
    string state_fname = folder_name + "state_ref_log.txt";
    string cntrl_fname = folder_name + "control_ref_log.txt";
    fstream state_fstrm(state_fname);
    fstream cntrl_fstrm(cntrl_fname);
    IOFormat eigFormat(Eigen::FullPrecision, 0, ",");

    if (!state_fstrm.is_open())
    {
        printf("Failed to open state reference log or file does not exist \n");
        return;
    }
    if (!cntrl_fstrm.is_open())
    {
        printf("Failed to open control reference log or file does not exist \n");
        return;
    }
    for (size_t i = 0; i < contactTraj.size(); i++)
    {
        int horizon = contactTraj[i].length() - 1;
        for (size_t k = 0; k < horizon; k++)
        {
            cntrl_fstrm << contactTraj[i].U[k].transpose().format(eigFormat) << std::endl;
            state_fstrm << contactTraj[i].X[k].transpose().format(eigFormat) << std::endl;
        }
        state_fstrm << contactTraj[i].X[horizon].transpose().format(eigFormat) << std::endl;
    }
}

template class ContactSchedule<ModelInfo<double, 24, 24, 0>>;