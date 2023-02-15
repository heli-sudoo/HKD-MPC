#include <iostream>
#include <cassert>
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
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::set_initial_condition(DVec<T> &x0_)
{
    X->at(0) = x0_;
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
void SinglePhase<T, xs, us, ys>::get_value_info_at_init(T &dV_, DVec<T> &G_, DMat<T> &H_)
{
    dV_ = dV->at(0);
    G_ = G->at(0);
    H_ = H->at(0);
}

template <typename T, size_t xs, size_t us, size_t ys>
DVec<T> SinglePhase<T, xs, us, ys>::get_terminal_state()
{
    return X->back();
}

template <typename T, size_t xs, size_t us, size_t ys>
DVec<T> SinglePhase<T, xs, us, ys>::get_terminal_state_nominal()
{
    return Xbar->back();
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

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::forward_sweep(T eps, HSDDP_OPTION &option, bool calc_partial)
{
    V->at(0) = 0;
    T Vprev = 0;
    actual_cost = 0;
    VecM<T, xs> delta_x;
    VecM<T, us> delta_u;
    VecM<T, ys> delta_y;
    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    int k = 0;
    for (k = 0; k < phase_horizon; k++)
    {
        /* update control */
        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * (X->at(k) - Xbar->at(k));
        /* run dynamics */
        dynamics(X->at(k + 1), Y->at(k), X->at(k), U->at(k));
        if (calc_partial)
        { // This flag is turned off when performing line search for speed up
            dynamics_partial(A->at(k), B->at(k), C->at(k), D->at(k), X->at(k), U->at(k));
        }
        /* compute running cost*/
        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);
        if (calc_partial)
        {
            /* compute running cost*/
            costContainer.running_cost_par(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);
        }

        /* compute path constraints */
        constraintContainer.compute_path_constraints(X->at(k), U->at(k), Y->at(k), k);
        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_with_pconstr(rcostData->at(k), pconstrsData, reb_params, calc_partial);
        }
        // V->at(k) = Vprev + rcostData->at(k).l;
        // Vprev = V->at(k);
        actual_cost += rcostData->at(k).l;
    }
    /* compute terminal cost and its partials */
    costContainer.terminal_cost(*tcostData, X->at(k), k);
    if (calc_partial)
    {
        // costContainer.terminal_cost_par(*tcostData, delta_x);
        costContainer.terminal_cost_par(*tcostData, X->at(k), k);
    }
    /* compute terminal constraint */
    constraintContainer.compute_terminal_constraints(X->at(k));
    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_with_tconstr(tconstrsData, al_params, calc_partial);
    }
    // V->at(phase_horizon) = Vprev + tcostData->Phi;
    actual_cost += tcostData->Phi;
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

    T exp_cost_change_1(0);
    T exp_cost_change_2(0);
    
    dX->at(0) = dx_init;
    for (int k = 0; k < phase_horizon; k++)
    {
        auto &qk = rcostData->at(k).lx;
        auto &rk = rcostData->at(k).lu;
        auto &Qk = rcostData->at(k).lxx;
        auto &Rk = rcostData->at(k).luu;
        auto &Pk = rcostData->at(k).lux;
        auto &Ak = A->at(k);
        auto &Bk = B->at(k);
        auto &dxk = dX->at(k);
        auto &Kk = K->at(k);
        auto &uff_k = dU->at(k);

        duk.setZero();
        duk = eps * uff_k + Kk * dxk;
        dX->at(k + 1) = Ak * dxk + Bk * duk + eps * Defect->at(k + 1);

        exp_cost_change_1 += qk.dot(dxk) + rk.dot(duk);
        exp_cost_change_2 += 0.5 * dxk.transpose() * Qk * dxk;
        exp_cost_change_2 += 0.5 * duk.transpose() * Rk * duk;
        exp_cost_change_2 += 0.5 * duk.transpose() * Pk * dxk;
    }

    // Terminal state
    auto &dxk = dX->back();
    exp_cost_change_1 += tcostData->Phix.dot(dxk);
    exp_cost_change_2 += 0.5 * dxk.transpose() * tcostData->Phixx * dxk;
}

/* Equivalent to system roll-out when eps = 0 */
template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::hybrid_rollout(T eps, HSDDP_OPTION &option)
{
    VecM<T, xs> dxk;
    VecM<T, xs> xk_next;
    const auto &SS_set = option.SS_set; // set of the shooting state time index (relative to the current phase)

    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    // V->at(0) = 0;
    // T Vprev = 0;
    actual_cost = 0;
    X->at(0) = x_init;
    Xsim->at(0) = xsim_init;

    int k = 0;
    for (int k = 0; k < phase_horizon; k++)
    {
        dxk = X->at(k) - Xbar->at(k);

        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * dxk;

        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        dynamics(Xsim->at(k + 1), Y->at(k), X->at(k), U->at(k));                

        auto it = std::find(SS_set.begin(), SS_set.end(), k + 1);
        if (it == SS_set.end()) // if k is a roll-out state        
        {
            X->at(k + 1) = Xsim->at(k + 1);
        }
        else // if k is a shooting state
        {
            X->at(k + 1) = Xbar->at(k+1) + eps * dX->at(k+1);
        }
        
        /* compute running cost*/
        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        /* compute path constraints */
        constraintContainer.compute_path_constraints(X->at(k), U->at(k), Y->at(k), k);

        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_with_pconstr(rcostData->at(k), pconstrsData, reb_params, 0);
        }

        // V->at(k) = Vprev + rcostData->at(k).l;
        // Vprev = V->at(k);
        actual_cost += rcostData->at(k).l;
    }

     /* compute terminal cost and its partials */
    costContainer.terminal_cost(*tcostData, X->at(k), k);
    
    /* compute terminal constraint */
    constraintContainer.compute_terminal_constraints(X->at(k));

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_with_tconstr(tconstrsData, al_params, 0);
    }
    // V->at(phase_horizon) = Vprev + tcostData->Phi;
    actual_cost = tcostData->Phi;    
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::nonlinear_rollout(T eps, HSDDP_OPTION &option)
{
    VecM<T, xs> dxk;
    VecM<T, xs> xk_next;
    const auto &SS_set = option.SS_set; // set of the shooting state time index (relative to the current phase)

    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    // V->at(0) = 0;
    // T Vprev = 0;

    actual_cost = 0;
    X->at(0) = x_init;
    Xsim->at(0) = xsim_init;

    int k = 0;
    for (int k = 0; k < phase_horizon; k++)
    {
        dxk = X->at(k) - Xbar->at(k);

        U->at(k) = Ubar->at(k) + eps * dU->at(k) + K->at(k) * dxk;

        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        dynamics(Xsim->at(k + 1), Y->at(k), X->at(k), U->at(k));        
        
        X->at(k + 1) = Xsim->at(k + 1) - (1 - eps) * Defect->at(k+1);

        /* compute running cost*/
        costContainer.running_cost(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        /* compute path constraints */
        constraintContainer.compute_path_constraints(X->at(k), U->at(k), Y->at(k), k);

        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_with_pconstr(rcostData->at(k), pconstrsData, reb_params, 0);
        }

        // V->at(k) = Vprev + rcostData->at(k).l;
        // Vprev = V->at(k);
        actual_cost += rcostData->at(k).l;
    }

     /* compute terminal cost and its partials */
    costContainer.terminal_cost(*tcostData, X->at(k), k);
    
    /* compute terminal constraint */
    constraintContainer.compute_terminal_constraints(X->at(k));

    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_with_tconstr(tconstrsData, al_params, 0);
    }
    // V->at(phase_horizon) = Vprev + tcostData->Phi;
    actual_cost += tcostData->Phi;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::LQ_approximation(HSDDP_OPTION& option)
{
    vector<IneqConstrData<T, xs, us, ys>> pconstrsData;
    vector<TConstrData<T, xs>> tconstrsData;
    vector<REB_Param_Struct<T>> reb_params;
    vector<AL_Param_Struct<T>> al_params;

    /* LQ approximation for all intermediate states */
    for (int k = 0; k < phase_horizon; k++)
    {
        dynamics_partial(A->at(k), B->at(k), C->at(k), D->at(k), X->at(k), U->at(k));

        costContainer.running_cost_par(rcostData->at(k), X->at(k), U->at(k), Y->at(k), dt, k);

        /* update running cost with path constraints using ReB method */
        if (option.ReB_active)
        {
            constraintContainer.get_path_constraints(pconstrsData, k);
            constraintContainer.get_reb_params(reb_params, k);
            update_running_cost_with_pconstr(rcostData->at(k), pconstrsData, reb_params, 1);
        }
    }
    
    /* update terminal cost with terminal constraint using AL */
    if (option.AL_active)
    {
        constraintContainer.get_terminal_constraints(tconstrsData);
        constraintContainer.get_al_params(al_params);
        update_terminal_cost_with_tconstr(tconstrsData, al_params, 1);
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
bool SinglePhase<T, xs, us, ys>::backward_sweep(T regularization, T dVprime, DVec<T> Gprime, DMat<T> Hprime)
{
    bool success = true;
    dV->at(phase_horizon) = dVprime;
    G->at(phase_horizon) = tcostData->Phix + Gprime;
    H->at(phase_horizon) = tcostData->Phixx + Hprime;

    for (int k = phase_horizon - 1; k >= 0; k--)
    {
        // compute Q function approximation
        auto &rcost = rcostData->at(k);
        auto &Ak = A->at(k);
        auto &Bk = B->at(k);
        auto &Ck = C->at(k);
        auto &Dk = D->at(k);
        auto &Gnext = G->at(k + 1);
        auto &Hnext = H->at(k + 1);

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
        dV->at(k) = dV->at(k + 1) - Qu.transpose() * Quu.inverse() * Qu;
    }
    return success;
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_pconstr(RCostData<T, xs, us, ys> &rcost,
                                                                  vector<IneqConstrData<T, xs, us, ys>> &pconstrsData,
                                                                  vector<REB_Param_Struct<T>> &reb_params,
                                                                  bool flag)
{
    compute_barrier(pconstrsData, reb_params); // update barrier data B, Bd, Bdd
    for (size_t i = 0; i < pConstrs_size; i++)
    {
        const auto &eps = reb_params[i].eps;
        const auto &c = pconstrsData[i];

        rcost.l += eps * Bar[i] * dt;
        if (flag)
        {
            rcost.lx += eps * Bard[i] * c.gx * dt;
            rcost.lu += eps * Bard[i] * c.gu * dt;
            rcost.ly += eps * Bard[i] * c.gy * dt;
            rcost.lxx += eps * dt * (Bardd[i] * c.gx * c.gx.transpose() + Bard[i] * c.gxx);
            rcost.luu += eps * dt * (Bardd[i] * c.gu * c.gu.transpose() + Bard[i] * c.guu);
            rcost.lyy += eps * dt * (Bardd[i] * c.gy * c.gy.transpose() + Bard[i] * c.gyy);
        }
    }
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_running_cost_with_smooth()
{
}

template <typename T, size_t xs, size_t us, size_t ys>
void SinglePhase<T, xs, us, ys>::update_terminal_cost_with_tconstr(vector<TConstrData<T, xs>> &tconstrsData,
                                                                   vector<AL_Param_Struct<T>> &al_params,
                                                                   bool flag)
{
    for (size_t i = 0; i < tConstrs_size; i++)
    {
        const T &sigma = al_params[i].sigma;
        const T &lambda = al_params[i].lambda;
        const auto &e = tconstrsData[i];
        tcostData->Phi += .5 * sigma * pow(e.h, 2) + lambda * e.h;
        if (flag)
        {
            tcostData->Phix += sigma * e.hx * e.h + lambda * e.hx;
            // Gauss Newton method to approximate the exact hessian hxx = hx*hx.transpose
            tcostData->Phixx += sigma * (e.hx * e.hx.transpose() + e.h * e.hx * e.hx.transpose()) + lambda * e.hx * e.hx.transpose();
        }
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
    for (size_t i = 0; i < pConstrs_size; i++)
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
    traj->push_back_zero();
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

template class SinglePhase<double, 24, 24, 0>;