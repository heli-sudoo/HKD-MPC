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
    dV_1_out = dV_2;
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
    dX->at(0) = dx_init;
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
        dV_2 += 0.5 * dxk.transpose() * Qk * dxk;
        dV_2 += 0.5 * duk.transpose() * Rk * duk;
        dV_2 += 0.5 * duk.transpose() * Pk * dxk;
    }

    // Terminal state
    const auto &dxk = dX->back();
    dV_1 += tcostData->Phix.dot(dxk);
    dV_2 += 0.5 * dxk.transpose() * tcostData->Phixx * dxk;
}

/* Equivalent to system roll-out when eps = 0 */
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::hybrid_rollout(T eps, HSDDP_OPTION &option, bool is_last_phase)
{            
    int k = 0;
    X->at(0) = x_init;
    Xsim->at(0) = xsim_init;
    
    for (k = 0; k < phase_horizon; k++)
    {
        /* update control */
        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * (X->at(k) - Xbar->at(k));

        /* run dynamics */
        dynamics(Xsim->at(k + 1), Y->at(k), X->at(k), U->at(k));

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
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::nonlinear_rollout(T eps, HSDDP_OPTION &option)
{
    VecM<T, xs> dxk;
    VecM<T, xs> xk_next;    

    int k = 0;
    actual_cost = 0;
    X->at(0) = x_init;
    Xsim->at(0) = xsim_init;    
    for (int k = 0; k < phase_horizon; k++)
    {
        dxk = X->at(k) - Xbar->at(k);

        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * dxk;

        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        dynamics(Xsim->at(k + 1), Y->at(k), X->at(k), U->at(k));

        X->at(k + 1) = Xsim->at(k + 1) - (1 - eps) * Defect_bar->at(k + 1);

    }

    /* compute terminal constraint */
    constraintContainer.compute_terminal_constraints(X->at(k));

    /* compute defects */    
    compute_defect();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::compute_cost(const HSDDP_OPTION& option)
{    
    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    int k =0;
    actual_cost = 0;

    for (k = 0; k < phase_horizon; k++)
    {
        /* compute running cost*/
        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);
        
        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_with_pconstr(rcostData->at(k), pconstrsData, reb_params);
        }
        actual_cost += rcostData->at(k).l;           
    }
    
    /* compute terminal cost */
    costContainer.terminal_cost(*tcostData, X->at(k), k);

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_with_tconstr(tconstrsData, al_params);
    }
    actual_cost += tcostData->Phi;   
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::LQ_approximation(HSDDP_OPTION &option)
{
    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    int k = 0;
    /* LQ approximation for all intermediate states */
    for (k = 0; k < phase_horizon; k++)
    {   
        /* compute dynamics partial*/
        dynamics_partial(A->at(k), B->at(k), C->at(k), D->at(k), X->at(k), U->at(k));

        /* compute running cost partial*/
        costContainer.running_cost_par(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        /* compute path constraints partial*/
        constraintContainer.compute_path_constraints_par(X->at(k), U->at(k), Y->at(k), k);

        /* update running cost partial with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_par_with_pconstr(rcostData->at(k), pconstrsData, reb_params);
        }
    }

    /* compute terminal cost partial*/
    costContainer.terminal_cost_par(*tcostData, X->at(k), k);

    /* compute terminal constraint partial*/
    constraintContainer.compute_terminal_constraints_par(X->at(k));

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_par_with_tconstr(tconstrsData, al_params);
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
        auto &Gnext = G->at(k + 1);

        // account for the defect
        Gnext += Hnext * Defect->at(k + 1);

        // standard update equations for Q
        Qx = rcost.lx + Ak.transpose() * Gnext + Ck.transpose() * rcost.ly;
        Qu = rcost.lu + Bk.transpose() * Gnext + Dk.transpose() * rcost.ly;
        Qxx = rcost.lxx + Ck.transpose() * rcost.lyy * Ck + Ak.transpose() * Hnext * Ak;
        Quu = rcost.luu + Dk.transpose() * rcost.lyy * Dk + Bk.transpose() * Hnext * Bk;
        Qux = rcost.lux + Dk.transpose() * rcost.lyy * Ck + Bk.transpose() * Hnext * Ak;

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
    return success;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_pconstr(RCostData<T, xs, us, ys> &rcost,
                                                                  vector<IneqConstrData<T, xs, us, ys>> &pconstrsData,
                                                                  vector<REB_Param_Struct<T>> &reb_params,
                                                                  int flag)
{
    if (flag == 1)
    {
        compute_barrier(pconstrsData, reb_params); // update barrier data B, Bd, Bdd
    }

    for (int i = 0; i < pConstrs_size; i++)
    {
        const auto &eps = reb_params[i].eps;
        const auto &c = pconstrsData[i];

        switch (flag)
        {
        case 1:
            rcost.l += eps * Bar[i] * dt;
            break;
        case 2:
            rcost.lx += eps * Bard[i] * c.gx * dt;
            rcost.lu += eps * Bard[i] * c.gu * dt;
            rcost.ly += eps * Bard[i] * c.gy * dt;
            rcost.lxx += eps * dt * (Bardd[i] * c.gx * c.gx.transpose() + Bard[i] * c.gxx);
            rcost.luu += eps * dt * (Bardd[i] * c.gu * c.gu.transpose() + Bard[i] * c.guu);
            rcost.lyy += eps * dt * (Bardd[i] * c.gy * c.gy.transpose() + Bard[i] * c.gyy);
            break;
        case 3:
            rcost.l += eps * Bar[i] * dt;
            rcost.lx += eps * Bard[i] * c.gx * dt;
            rcost.lu += eps * Bard[i] * c.gu * dt;
            rcost.ly += eps * Bard[i] * c.gy * dt;
            rcost.lxx += eps * dt * (Bardd[i] * c.gx * c.gx.transpose() + Bard[i] * c.gxx);
            rcost.luu += eps * dt * (Bardd[i] * c.gu * c.gu.transpose() + Bard[i] * c.guu);
            rcost.lyy += eps * dt * (Bardd[i] * c.gy * c.gy.transpose() + Bard[i] * c.gyy);
            break;
        }
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_pconstr(RCostData<T, xs, us, ys> &rcost,
                                                                  vector<IneqConstrData<T, xs, us, ys>> &pconstrsData,
                                                                  vector<REB_Param_Struct<T>> &reb_params)
{
    compute_barrier(pconstrsData, reb_params); // update barrier data B, Bd, Bdd

    for (int i = 0; i < pConstrs_size; i++)
    {
        const auto &eps = reb_params[i].eps;
        const auto &c = pconstrsData[i];
        rcost.l += eps * Bar[i] * dt;
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_par_with_pconstr(RCostData<T, xs, us, ys> &rcost,
                                                                      vector<IneqConstrData<T, xs, us, ys>> &pconstrsData,
                                                                      vector<REB_Param_Struct<T>> &reb_params)
{

    for (int i = 0; i < pConstrs_size; i++)
    {
        const auto &eps = reb_params[i].eps;
        const auto &c = pconstrsData[i];
        rcost.lx += eps * Bard[i] * c.gx * dt;
        rcost.lu += eps * Bard[i] * c.gu * dt;
        rcost.ly += eps * Bard[i] * c.gy * dt;
        rcost.lxx += eps * dt * (Bardd[i] * c.gx * c.gx.transpose() + Bard[i] * c.gxx);
        rcost.luu += eps * dt * (Bardd[i] * c.gu * c.gu.transpose() + Bard[i] * c.guu);
        rcost.lyy += eps * dt * (Bardd[i] * c.gy * c.gy.transpose() + Bard[i] * c.gyy);
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_smooth()
{
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_with_tconstr(vector<TConstrData<T, xs>> &tconstrsData,
                                                                   vector<AL_Param_Struct<T>> &al_params,
                                                                   int flag)
{
    for (int i = 0; i < tConstrs_size; i++)
    {
        const T &sigma = al_params[i].sigma;
        const T &lambda = al_params[i].lambda;
        const auto &e = tconstrsData[i];

        switch (flag)
        {
        case 1:
            tcostData->Phi += .5 * sigma * pow(e.h, 2) + lambda * e.h;
            break;
        case 2:
            tcostData->Phix += sigma * e.hx * e.h + lambda * e.hx;
            // Gauss Newton method to approximate the exact hessian hxx = hx*hx.transpose
            tcostData->Phixx += sigma * (e.hx * e.hx.transpose() + e.h * e.hx * e.hx.transpose()) + lambda * e.hx * e.hx.transpose();
            break;
        case 3:
            tcostData->Phi += .5 * sigma * pow(e.h, 2) + lambda * e.h;
            tcostData->Phix += sigma * e.hx * e.h + lambda * e.hx;
            // Gauss Newton method to approximate the exact hessian hxx = hx*hx.transpose
            tcostData->Phixx += sigma * (e.hx * e.hx.transpose() + e.h * e.hx * e.hx.transpose()) + lambda * e.hx * e.hx.transpose();
            break;
        }
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_with_tconstr(vector<TConstrData<T, xs>> &tconstrsData,
                                                                   vector<AL_Param_Struct<T>> &al_params)
{
    for (int i = 0; i < tConstrs_size; i++)
    {
        const T &sigma = al_params[i].sigma;
        const T &lambda = al_params[i].lambda;
        const auto &e = tconstrsData[i];

        tcostData->Phi += .5 * sigma * pow(e.h, 2) + lambda * e.h;
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_par_with_tconstr(vector<TConstrData<T, xs>> &tconstrsData,
                                                                       vector<AL_Param_Struct<T>> &al_params)
{
    for (int i = 0; i < tConstrs_size; i++)
    {
        const T &sigma = al_params[i].sigma;
        const T &lambda = al_params[i].lambda;
        const auto &e = tconstrsData[i];

        tcostData->Phix += sigma * e.hx * e.h + lambda * e.hx;
        // Gauss Newton method to approximate the exact hessian hxx = hx*hx.transpose
        tcostData->Phixx += sigma * (e.hx * e.hx.transpose() + e.h * e.hx * e.hx.transpose()) + lambda * e.hx * e.hx.transpose();
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_nominal_trajectory()
{
    traj->update_nominal_vals();
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::compute_barrier(vector<IneqConstrData<T, xs, us, ys>> &pconstrsData,
                                                 vector<REB_Param_Struct<T>> &reb_params)
{
    Bar.setZero();
    Bard.setZero();
    Bardd.setZero();

    int k = 2; // order of approximating polynomial
    for (int i = 0; i < pConstrs_size; i++)
    {
        const auto &c = pconstrsData[i];
        const auto &delta = reb_params[i].delta;

        if (c.g > delta)
        {
            Bar[i] = -log(c.g);
            Bard[i] = -1.0 / c.g;
            Bardd[i] = pow(c.g, -2);
        }
        else
        {
            Bar[i] = .5 * (((c.g - 2 * delta) / delta) * ((c.g - 2 * delta) / delta) - 1) - log(delta);
            Bard[i] = (c.g - 2 * delta) / delta / delta;
            Bardd[i] = pow(delta, -2);
        }
    }
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
void SinglePhase<T, xs, us, ys>::push_back()
{
    // traj->push_back_zero();
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

template class SinglePhase<double, 24, 24, 0>;