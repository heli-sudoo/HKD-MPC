#ifndef HKDLOG_H
#define HKDLOG_H
#include <fstream>
#include "HSDDP_CPPTypes.h"
#include "HKDProblem.h"
#include "TrajectoryManagement.h"
#include <iomanip>      // std::setprecision

using Eigen::IOFormat;
template <class TrajPtrs>
void log_trajectory_sequence(TrajPtrs &traj_seq)
{
    string folder_name = "../log/";
    string state_log_fname = folder_name + "state_log.txt";
    string cntrl_log_fname = folder_name + "control_log.txt";
    string cost_log_fname = folder_name + "cost_log.txt";
    string value_grad_fname = folder_name + "value_grad_log.txt";

    fstream state_fstrm(state_log_fname, std::ios_base::out | std::ios_base::trunc);
    fstream cntrl_fstrm(cntrl_log_fname, std::ios_base::out | std::ios_base::trunc);
    fstream cost_fstrm(cost_log_fname, std::ios_base::out | std::ios_base::trunc);
    fstream value_grad_fstrm(value_grad_fname, std::ios_base::out | std::ios_base::trunc);

    if (!state_fstrm.is_open())
    {
        printf("Failed to open state log file or file not exists \n");
        return;
    }
    if (!cntrl_fstrm.is_open())
    {
        printf("Failed to open control log file or file not exists \n");
        return;
    }
    if (!cost_fstrm.is_open())
    {
        printf("Failed to open cost log file or file not exists \n");
        return;
    }
    if (!value_grad_fstrm.is_open())
    {
        printf("Failed to open value_grad log file or file not exists \n");
        return;
    }

    state_fstrm.clear();
    cntrl_fstrm.clear();
    cost_fstrm.clear();
    value_grad_fstrm.clear();
    printf("Start logging trajectory \n");
    for (auto traj:traj_seq)
    {
        size_t horizon = traj->horizon;
        size_t k(0);
        for (k = 0; k < horizon; k++)
        {
            cntrl_fstrm << eigenToString(traj->Ubar[k].transpose()) << std::endl;
            state_fstrm << eigenToString(traj->Xbar[k].transpose()) << std::endl;
            value_grad_fstrm << eigenToString(traj->G[k].transpose()) << std::endl;
            cost_fstrm << traj->rcostData[k].l << std::endl;
        }
        cntrl_fstrm << eigenToString(traj->Ubar[k].transpose()) << std::endl;
        state_fstrm << eigenToString(traj->Xbar[k].transpose()) << std::endl; // log terminal state
        value_grad_fstrm << eigenToString(traj->G[k].transpose()) << std::endl;
        cost_fstrm  << traj->tcostData.Phi << std::endl;
    }
    printf("DONE! \n");
    state_fstrm.close();
    cntrl_fstrm.close();
    cost_fstrm.close();
    value_grad_fstrm.close();
}

template <typename T>
void log_contact_sequence(deque<VecM<int, 4>>& ctactSeq, deque<T>& horizons)
{
    string folder_name = "../log/";
    string contact_log_fname = folder_name + "contact_log.txt";

    fstream contact_fstrm(contact_log_fname);
    if (!contact_fstrm.is_open())
    {
        printf("Failed to open contact log file or file not exists \n");
        return;
    }
    contact_fstrm.clear();
    printf("Start logging contact trajectory \n");
    for (int i = 0; i < horizons.size(); i++)
    {
        for (int k = 0; k <= horizons[i]; k++)
        {
            contact_fstrm << ctactSeq[i].transpose() << std::endl;
        }
        
    }    
    printf("DONE! \n");
    contact_fstrm.close();
}

#endif