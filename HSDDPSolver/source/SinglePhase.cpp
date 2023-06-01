#include <iostream>
#include <cassert>
#include <numeric>
#include "SinglePhase.h"
#include "HSDDP_Utils.h"

template <typename T, size_t xs, size_t us, size_t ys>
SinglePhase<T, xs, us, ys>::SinglePhase()
{
    Qx.setZero();
    Qu.setZero();
    Qxx.setZero();
    Quu.setZero();
    Qux.setZero();
    Ixx.setIdentity();
    Iuu.setIdentity();
    Quu_inv.setZero();
    tConstrs_size = 0;
    pConstrs_size = 0;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::initialization()
{
    pConstrs_size = constraintContainer.pcontrs_size;
    tConstrs_size = constraintContainer.tconstr_size;
    Bar.setZero(pConstrs_size);
    Bard.setZero(pConstrs_size);
    Bardd.setZero(pConstrs_size);

    x_init.setZero();
    xsim_init.setZero();
    dx_init.setZero();    
    SS_set.clear();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_initial_condition(DVec<T> &x0_)
{
    x_init = x0_;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_initial_condition(DVec<T> &x0_, DVec<T> &x_sim_0_)
{
    x_init = x0_;
    xsim_init = x_sim_0_;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_initial_condition_dx(DVec<T> &dx0_)
{
    dx_init = dx0_;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_nominal_initial_condition(DVec<T> &x0_)
{
    Xbar->at(0) = x0_;
}
template <typename T, size_t xs, size_t us, size_t ys>
DVec<T> SinglePhase<T, xs, us, ys>::resetmap(DVec<T> &x_)
{
    DVec<T> xnext;
    if (resetmap_func_handle != nullptr)
    {
        resetmap_func_handle(xnext, x_);
    }
    return xnext;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::resetmap_partial(DMat<T> &Px_, DVec<T> &x_)
{
    if (resetmap_partial_func_handle != nullptr)
    {
        resetmap_partial_func_handle(Px_, x_);
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_value_approx(DVec<T> &G_out, DMat<T> &H_out)
{
    G_out = G->at(0);
    H_out = H->at(0);
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_exp_cost_change(T &dV_1_out, T &dV_2_out)
{
    dV_1_out = dV_1;
    dV_2_out = dV_2;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_terminal_state(DVec<T> &xend_out)
{
    xend_out = X->back();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_terminal_state(DVec<T> &xend_out, DVec<T> &xsim_end_out)
{
    xend_out = X->back();
    xsim_end_out = Xsim->back();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_terminal_state_dx(DVec<T> &dx_end_out)
{
    dx_end_out = dX->back();
}

template <typename T, size_t xs, size_t us, size_t ys>
T SinglePhase<T, xs, us, ys>::get_actual_cost()
{
    return actual_cost;
}

template <typename T, size_t xs, size_t us, size_t ys>
T SinglePhase<T, xs, us, ys>::get_max_pconstrs()
{
    return constraintContainer.get_max_pconstrs();
}

template <typename T, size_t xs, size_t us, size_t ys>
T SinglePhase<T, xs, us, ys>::get_max_tconstrs()
{
    return constraintContainer.get_max_tconstrs();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_trajectory(shared_ptr<Trajectory<T, xs, us, ys>> traj_)
{
    traj = traj_;
    update_trajectory_ptrs();
}

/*
brief:
    Obtain the search direction dX for the shooting state using the linearized dynaamics
    Compute the expected cost change
*/
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::linear_rollout(T eps, HSDDP_OPTION &option)
{
    VecM<T, us> duk;

    dV_1 = 0;
    dV_2 = 0;
    dX->at(0) = dx_init + eps * Defect->at(0);
    for (int k = 0; k < phase_horizon; k++)
    {
        const auto &qk = rcostData->at(k).lx;
        const auto &rk = rcostData->at(k).lu;
        const auto &Qk = rcostData->at(k).lxx;
        const auto &Rk = rcostData->at(k).luu;
        const auto &Pk = rcostData->at(k).lux;
        const auto &Ak = A->at(k);
        const auto &Bk = B->at(k);
        const auto &dxk = dX->at(k);
        const auto &Kk = K->at(k);
        const auto &uff_k = dU->at(k);

        duk = eps * uff_k + Kk * dxk;
        dX->at(k + 1) = Ak * dxk + Bk * duk + eps * Defect->at(k + 1);

        dV_1 += qk.dot(dxk) + rk.dot(duk);
        dV_2 += dxk.transpose() * Qk * dxk;
        dV_2 += duk.transpose() * Rk * duk;
        dV_2 += duk.transpose() * Pk * dxk;
    }

    // Terminal state
    const auto &dxk = dX->back();
    dV_1 += tcostData->Phix.dot(dxk);
    dV_2 += dxk.transpose() * tcostData->Phixx * dxk;
}

/* Equivalent to system roll-out when eps = 0 */
template <typename T, size_t xs, size_t us, size_t ys>
bool SinglePhase<T, xs, us, ys>::hybrid_rollout(T eps, HSDDP_OPTION &option, bool is_last_phase)
{    
    (void) (is_last_phase);               
    Xsim->at(0) = x_init;
    
    if (!SS_set.empty() && SS_set.front()==0)
    {
        X->at(0) = Xbar->at(0) + eps * dX->at(0);
    }else
    {
        X->at(0) = x_init;
    }
    
    
    int k = 0;
    for (k = 0; k < phase_horizon; k++)
    {
        /* update control */
        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * (X->at(k) - Xbar->at(k));

        /* run dynamics */
        dynamics(Xsim->at(k + 1), Y->at(k), X->at(k), U->at(k), t_offset + k*dt);

        if (Xsim->at(k+1).template lpNorm<2>() > 1e6)
        {
            return false;
        }
        
        /* check wether k+1 is a shooting state */
        auto it = std::find(SS_set.begin(), SS_set.end(), k + 1);

        /* If k+1 is a shooting state, and multiple shooting is set to true*/
        if ((option.MS) && (SS_set.size()>0) && (it != SS_set.end()))
        {
            X->at(k + 1) = Xbar->at(k + 1) + eps * dX->at(k + 1);
        }else
        {
            X->at(k + 1) = Xsim->at(k + 1);
        }

        /* compute path constraints */
        constraintContainer.compute_path_constraints(X->at(k), U->at(k), Y->at(k), k);          
    }    

    /* compute terminal constraint */
    constraintContainer.compute_terminal_constraints(X->at(k));

    /* compute defect */
    compute_defect();

    return true;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::compute_cost(const HSDDP_OPTION& option)
{       
    size_t k =0;
    actual_cost = 0;    
    for (k = 0; k < phase_horizon; k++)
    {
        /* compute running cost*/
        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, t_offset + k*dt);
        
        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {            
            update_running_cost_with_pconstr(k);            
        }
        actual_cost += rcostData->at(k).l;           
    }
        
    /* compute terminal cost */
    costContainer.terminal_cost(*tcostData, X->at(k), t_offset + k*dt);

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {        
        update_terminal_cost_with_tconstr();
    }
    actual_cost += tcostData->Phi;   
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::LQ_approximation(HSDDP_OPTION &option)
{    
    int k = 0;
    /* LQ approximation for all intermediate states */
    for (k = 0; k < phase_horizon; k++)
    {   
        /* compute dynamics partial*/
        dynamics_partial(A->at(k), B->at(k), C->at(k), D->at(k), X->at(k), U->at(k), t_offset + k*dt);

        /* compute running cost partial*/
        costContainer.running_cost_par(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, t_offset + k*dt);

        /* update running cost partial with path constraints using ReB method */
        if (option.ReB_active)
        {
             /* compute path constraints partial*/
            constraintContainer.compute_path_constraints_par(X->at(k), U->at(k), Y->at(k), k);    
            update_running_cost_par_with_pconstr(k);
        }
    }

    /* compute terminal cost partial*/
    costContainer.terminal_cost_par(*tcostData, X->at(k), t_offset + k*dt);

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
         /* compute terminal constraint partial*/
        constraintContainer.compute_terminal_constraints_par(X->at(k));                
        update_terminal_cost_par_with_tconstr();
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
bool SinglePhase<T, xs, us, ys>::backward_sweep(T regularization, DVec<T> Gprime, DMat<T> Hprime)
{
    bool success = true;
    G->at(phase_horizon) = tcostData->Phix + Gprime;
    H->at(phase_horizon) = tcostData->Phixx + Hprime;

    dV_1 = 0;
    dV_2 = 0;
    for (int k = phase_horizon - 1; k >= 0; k--)
    {
        // compute Q function approximation
        const auto &rcost = rcostData->at(k);
        const auto &Ak = A->at(k);
        const auto &Bk = B->at(k);
        const auto &Ck = C->at(k);
        const auto &Dk = D->at(k);
        const auto &Hnext = H->at(k + 1);
        const auto &Defect_next = Defect->at(k + 1);
        auto Gnext = G->at(k + 1);

        // account for the defect
        Gnext += Hnext * Defect_next;

        // standard update equations for Q
        Qx = rcost.lx + Ak.transpose() * Gnext;
        Qu = rcost.lu + Bk.transpose() * Gnext;
        Qxx = rcost.lxx + Ak.transpose() * Hnext * Ak;
        Quu = rcost.luu + Bk.transpose() * Hnext * Bk;
        Qux = rcost.lux + Bk.transpose() * Hnext * Ak;

        if (ys > 0)
        {
            Qx += Ck.transpose() * rcost.ly;
            Qu += Dk.transpose() * rcost.ly;
            Qxx += Ck.transpose() * rcost.lyy * Ck;
            Quu += Dk.transpose() * rcost.lyy * Dk;
            Qux += Dk.transpose() * rcost.lyy * Ck;
        }
        

        // regularizatoin
        Qxx += Ixx * regularization;
        Quu += Iuu * regularization;
        Quu_chol = Quu_chol.compute(Quu - Iuu * 1e-9); // Cholesky decomposition of Quu
        // If Quu not PSD, break and return false
        if (!Quu_chol.isPositive())
        {
            success = false;
            break;
        }
        // compute value function approximation and sum up expected cost change at each time step
        // Symmetrize Quu_inv and Qxx. Numerical issue would occur otherwise.
        Quu_inv = (Quu.inverse() + Quu.inverse().transpose()) / 2;
        Qxx = (Qxx + Qxx.transpose()) / 2;

        // compute value function approximation
        dU->at(k) = -Quu_inv * Qu;
        K->at(k) = -Quu_inv * Qux;
        G->at(k) = Qx - Qux.transpose() * Quu_inv * Qu;
        H->at(k) = Qxx - Qux.transpose() * Quu_inv * Qux;
        T dV_k = -Qu.transpose() * dU->at(k);

        dV_1 -= dV_k;
        dV_2 += dV_k;
    }   
    // Update the value gradient using the defect at the initial condition
    G->at(0) += H->at(0) * Defect->at(0);
    return success;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_pconstr(size_t k)
{
    const auto& pConstraints = constraintContainer.pathConstraints;
    for (size_t i = 0; i < pConstraints.size(); i++)
    {
        pConstraints[i]->compute_ReB_cost(k);
        rcostData->at(k).l += dt * pConstraints[i]->ReB_cost;
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_par_with_pconstr(size_t k)
{
    const auto& pConstraints = constraintContainer.pathConstraints;
    for (size_t i = 0; i < pConstraints.size(); i++)
    {
        pConstraints[i]->compute_ReB_partials(k);
        rcostData->at(k).lu += dt * pConstraints[i]->ReB_grad_u;
        rcostData->at(k).lx += dt * pConstraints[i]->ReB_grad_x;
        rcostData->at(k).ly += dt * pConstraints[i]->ReB_grad_y;
        rcostData->at(k).luu += dt * pConstraints[i]->ReB_hess_u;
        rcostData->at(k).lxx += dt * pConstraints[i]->ReB_hess_x;
        rcostData->at(k).lyy += dt * pConstraints[i]->ReB_hess_y;
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_smooth()
{
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_with_tconstr()
{
    const auto& tConstraints = constraintContainer.terminalConstraints;
    for (size_t i = 0; i < tConstraints.size(); i++)
    {
        tConstraints[i]->compute_AL_cost();
        const T& AL_cost = tConstraints[i]->AL_cost;
        tcostData->Phi += AL_cost;
    }    
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_par_with_tconstr()
{
    const auto& tConstraints = constraintContainer.terminalConstraints;
    for (size_t i = 0; i < tConstraints.size(); i++)
    {
        tConstraints[i]->compute_AL_partials();
        const auto& AL_gradient = tConstraints[i]->AL_gradient;
        const auto& AL_hessian = tConstraints[i]->AL_hessian;

        tcostData->Phix += AL_gradient;
        tcostData->Phixx += AL_hessian;
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_nominal_trajectory()
{
    traj->update_nominal_vals();
}


template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_trajectory_ptrs()
{
    if (traj == nullptr)
    {
        printf("Pointer to phase trajectory is null \n");
        return;
    }
    Xbar = &(traj->Xbar);
    X = &(traj->X);
    Ubar = &(traj->Ubar);
    U = &(traj->U);
    Y = &(traj->Y);

    A = &(traj->A);
    B = &(traj->B);
    C = &(traj->C);
    D = &(traj->D);

    Xsim = &(traj->Xsim);
    Defect = &(traj->Defect);
    Defect_bar = &(traj->Defect_bar);
    dX = &(traj->dX);

    V = &(traj->V);
    dV = &(traj->dV);
    dU = &(traj->dU);
    K = &(traj->K);
    G = &(traj->G);
    H = &(traj->H);
    rcostData = &(traj->rcostData);
    tcostData = &(traj->tcostData);
    dt = traj->timeStep;
    phase_horizon = traj->horizon;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_AL_params(HSDDP_OPTION &option)
{
    constraintContainer.update_al_params(option);
}
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_REB_params(HSDDP_OPTION &option)
{
    constraintContainer.update_reb_params(option);
}

/*
    @brief  Push back n (zero) elements to the trajectory and path constraints
*/
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::push_back_default()
{    
    traj->push_back_state(traj->X.back());    
    constraintContainer.push_back_n(1);
    phase_horizon = traj->horizon;
}
/*
    @brief Remove n elements in the front from the trajectory and path constraints
*/
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::pop_front()
{
    traj->pop_front();
    constraintContainer.pop_front_n(1);
    phase_horizon = traj->horizon;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::get_trajectory(std::vector< std::vector<float>> & x_tau, 
                                                std::vector< std::vector<float>> & u_tau)
{
    std::vector<float> xk;
    std::vector<float> uk;
    for (int k = 0; k < phase_horizon; k++)
    {
        xk.clear();
        uk.clear();
       for (int i = 0; i < xs; i++)
       {
           xk.push_back(X->at(k)[i]);
           uk.push_back(U->at(k)[i]);
       }
       x_tau.push_back(xk);
       u_tau.push_back(uk);
    }    
}      

/*
@brief: Print the phase information (primarily used for debugging)  
*/
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::print()
{
    printf("************* Phase Information ***************\n");
    printf("           phase_horizon = %d \n", phase_horizon);
    printf("                      dt = %f \n", dt);
    printf("                t_offset = %f \n", t_offset);
    printf("         trajectory size = %d \n", traj->size());
    printf("num terminal constraints = %lu \n", constraintContainer.num_terminal_constraints());
    printf("    num path constraints = %lu \n\n", constraintContainer.num_path_constraints());
}

template class SinglePhase<double, 24, 24, 0>;
template class SinglePhase<double, 12, 12, 0>;
template class SinglePhase<double, 36, 12, 12>;