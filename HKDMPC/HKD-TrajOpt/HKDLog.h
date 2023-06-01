#ifndef HKDLOG_H
#define HKDLOG_H
#include <fstream>
#include <iomanip>      // std::setprecision

template <class TrajPtrs>
void log_trajectory_sequence(const std::string& folder_name, TrajPtrs &traj_seq)
{    
    std::string state_log_fname = folder_name + "state_log.txt";
    std::string cntrl_log_fname = folder_name + "control_log.txt";
    std::string cost_log_fname = folder_name + "cost_log.txt";
    std::string value_grad_fname = folder_name + "value_grad_log.txt";

    std::fstream state_fstrm(state_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cntrl_fstrm(cntrl_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cost_fstrm(cost_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream value_grad_fstrm(value_grad_fname, std::ios_base::out | std::ios_base::trunc);

    Eigen::IOFormat eigFormat(Eigen::FullPrecision, 0, ",");
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
        for (size_t k = 0; k < horizon; k++)
        {
            cntrl_fstrm << traj->Ubar[k].transpose().format(eigFormat) << std::endl;
            state_fstrm << traj->Xbar[k].transpose().format(eigFormat) << std::endl;
            value_grad_fstrm << traj->G[k].transpose().format(eigFormat) << std::endl;
            cost_fstrm <<  traj->rcostData[k].l << std::endl;
        }
        cntrl_fstrm << traj->Ubar[horizon-1].transpose().format(eigFormat) << std::endl;
        state_fstrm << traj->Xbar[horizon].transpose().format(eigFormat) << std::endl; // log terminal state
        value_grad_fstrm << traj->G[horizon].transpose().format(eigFormat) << std::endl;
        cost_fstrm << traj->tcostData.Phi << std::endl;
    }
    printf("DONE! \n");
    state_fstrm.close();
    cntrl_fstrm.close();
    cost_fstrm.close();
    value_grad_fstrm.close();
}



#endif