#pragma once
#ifndef HSDDP_UTILS_H
#define HSDDP_UTILS_H
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

struct TIME_PER_ITERATION
{
    int DDP_iter = 0;
    int n_bws = 0;           // number of backward sweep iterations
    double time_bws = 0;     // backward sweep time
    int n_fit = 0;           // number of iterations in line search
    double time_fit = 0;     // line search time
    double time_partial = 0; // time for computing derivatives (dynamics and cost)
};

template <typename Derived, typename T>
std::vector<int> find_eigen(const Eigen::DenseBase<Derived> &v, const T &val)
{
    typedef typename Derived::Scalar val_type;
    std::vector<int> index;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == static_cast<val_type>(val))
        {
            index.push_back(i);
        }
    }
    return index;
}

template <typename SeqContainer>
void append_vector(SeqContainer &v1, const SeqContainer &v2)
{
    for (auto &val : v2)
    {
        v1.push_back(val);
    }
}

// Check if two floating-point number are almost equal
template <typename T1, typename T2>
bool approx_eq_scalar(T1 n1, T2 n2)
{
    float tol = 1e-6; // 1e-6 for float, 1e-8 for double
    float err = std::abs(n1 - n2);
    if (err <= tol)
    {
        return true;
    }
    return false;
}

// Check if n1 < n2 or n1 approx_eq n2
template <typename T1, typename T2>
bool approx_leq_scalar(T1 n1, T2 n2)
{
    if (n1 < n2 || approx_eq_scalar(n1, n2))
    {
        return true;
    }
    return false;
}

// Check if n1 > n2 or n1 approx_eq n2
template <typename T1, typename T2>
bool approx_geq_scalar(T1 n1, T2 n2)
{
    if (n1 > n2 || approx_eq_scalar(n1, n2))
    {
        return true;
    }
    return false;
}

template <class TrajPtrs>
void log_trajectory_sequence(const std::string &folder_name, const TrajPtrs &traj_seq)
{
    std::string state_log_fname = folder_name + "state_log.txt";
    std::string cntrl_log_fname = folder_name + "control_log.txt";
    std::string cost_log_fname = folder_name + "cost_log.txt";
    std::string value_grad_fname = folder_name + "value_grad_log.txt";

    std::fstream state_fstrm(state_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cntrl_fstrm(cntrl_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cost_fstrm(cost_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream value_grad_fstrm(value_grad_fname, std::ios_base::out | std::ios_base::trunc);

    Eigen::IOFormat eigFormat(5, 0, ",");
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

    for (const auto &traj : traj_seq)
    {
        size_t horizon = traj->horizon;
        for (size_t k = 0; k < horizon; k++)
        {
            cntrl_fstrm << traj->Ubar[k].transpose().format(eigFormat) << "\n";
            state_fstrm << traj->Xbar[k].transpose().format(eigFormat) << "\n";
            value_grad_fstrm << traj->G[k].transpose().format(eigFormat) << "\n";
            cost_fstrm << traj->rcostData[k].l << "\n";
        }
        cntrl_fstrm << traj->Ubar[horizon - 1].transpose().format(eigFormat) << "\n";
        state_fstrm << traj->Xbar[horizon].transpose().format(eigFormat) << "\n"; // log terminal state
        value_grad_fstrm << traj->G[horizon].transpose().format(eigFormat) << "\n";
        cost_fstrm << traj->tcostData.Phi << "\n";
    }

    state_fstrm.close();
    cntrl_fstrm.close();
    cost_fstrm.close();
    value_grad_fstrm.close();

    std::cout << "Logged trajectories to folder " << folder_name << "\n";
}

template <class TrajPtr>
void log_a_trajectory(const std::string &folder_name, const TrajPtr &traj)
{
    const std::string& state_log_fname = folder_name + "state_log.txt";
    const std::string& cntrl_log_fname = folder_name + "control_log.txt";
    const std::string& cost_log_fname = folder_name + "cost_log.txt";
    const std::string& value_grad_fname = folder_name + "value_grad_log.txt";
    const std::string& A_log_fname = folder_name + "dynamics_partial_A.txt";
    const std::string& B_log_fname = folder_name + "dynamics_partial_B.txt";


    std::fstream state_fstrm(state_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cntrl_fstrm(cntrl_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream cost_fstrm(cost_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream value_grad_fstrm(value_grad_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream A_fstrm(A_log_fname, std::ios_base::out | std::ios_base::trunc);
    std::fstream B_fstrm(B_log_fname, std::ios_base::out | std::ios_base::trunc);

    Eigen::IOFormat eigFormat(5, 1, ",");
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
    if (!B_fstrm.is_open())
    {
        printf("Failed to open dynamics_partial_B log file or file not exists \n");
        return;
    }
    if (!A_fstrm.is_open())
    {
        printf("Failed to open dynamics_partial_A log file or file not exists \n");
        return;
    }

    state_fstrm.clear();
    cntrl_fstrm.clear();
    cost_fstrm.clear();
    value_grad_fstrm.clear();
    A_fstrm.clear();
    B_fstrm.clear();
    
    size_t k = 0;
    for (k = 0; k < traj->size()-1; k++)
    {
        cntrl_fstrm << traj->Ubar[k].transpose().format(eigFormat) << "\n";
        state_fstrm << traj->Xbar[k].transpose().format(eigFormat) << "\n";
        value_grad_fstrm << traj->G[k].transpose().format(eigFormat) << "\n";
        cost_fstrm << traj->rcostData[k].l << "\n";
        A_fstrm << traj->A[k] << "\n\n";
        B_fstrm << traj->B[k] << "\n\n";
    }
    cntrl_fstrm << traj->Ubar[k].transpose().format(eigFormat) << "\n";
    state_fstrm << traj->Xbar[k].transpose().format(eigFormat) << "\n"; // log terminal state
    value_grad_fstrm << traj->G[k].transpose().format(eigFormat) << "\n";
    cost_fstrm << traj->tcostData.Phi << "\n";

    state_fstrm.close();
    cntrl_fstrm.close();
    cost_fstrm.close();
    value_grad_fstrm.close();

    std::cout << "Logged trajectories to folder " << folder_name << "\n";
}

std::ostream &print_time(std::ostream &os, std::vector<TIME_PER_ITERATION> &time_ddp);

#endif