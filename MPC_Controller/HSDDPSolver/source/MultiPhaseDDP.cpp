#include "MultiPhaseDDP.h"
#include <algorithm>
#include "HSDDP_CompoundTypes.h"
#include "cTypes.h"

#define DEBUG

#ifdef TIME_BENCHMARK
#include <chrono>
using namespace std::chrono;
using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;
static int fit_iter = 0;
#endif // TIME_BENCHMARK


/*
    @brief: perform linear rollout to compute search direction for shooting state
            use before running hybrid rollout
*/
template <typename T>
void MultiPhaseDDP<T>::linear_rollout(T eps, HSDDP_OPTION &option)
{
    DVec<T> dx_init, dx_end;
    DVec<T> xend, xsim_end;
    DMat<T> Px; // resetmap partial
    dx0.setZero();
    dx_init = dx0;

    for (int i = 0; i < n_phases; i++)
    {
        if (i > 0)
        {
            phases[i - 1]->get_terminal_state_dx(dx_end);
            phases[i - 1]->get_terminal_state(xend, xsim_end);
            phases[i - 1]->resetmap_partial(Px, xend);
            dx_init = Px * dx_end;
        }

        phases[i]->set_initial_condition_dx(dx_init);

        phases[i]->linear_rollout(eps, option);
    }
}

/*
    @brief: perform hybrid rollout
            nominal system rollout is obtained with eps = 0
*/
template <typename T>
void MultiPhaseDDP<T>::hybrid_rollout(T eps, HSDDP_OPTION &option)
{
    actual_cost = 0;
    max_pconstr = 0;
    max_tconstr = 0;
    DVec<T> xinit = x0; // initial condition for each phase
    DVec<T> xsim_init = x0;
    DVec<T> xend;
    DVec<T> xsim_end;
    bool is_last_phase = false;

    run_before_forward_sweep(); // currently not used   

    for (int i = 0; i < n_phases; i++)
    {
        if (i > 0)
        {
            // If not the first phase, run resetmap at the end of previous phase
            phases[i - 1]->get_terminal_state(xend, xsim_end);
            xinit = phases[i - 1]->resetmap(xend);
            xsim_init = phases[i - 1]->resetmap(xsim_end);
        }        

        phases[i]->set_initial_condition(xinit, xsim_init); // Set initial condition of current phase
        phases[i]->hybrid_rollout(eps, option, is_last_phase);             // run hybrid rollout for each phase

        // update the maximum constraint violations
        max_pconstr = std::min(max_pconstr, phases[i]->get_max_pconstrs()); // should have non-positive value
        max_tconstr = std::max(max_tconstr, phases[i]->get_max_tconstrs()); // should have non-negative value
    }    
}

template <typename T>
void MultiPhaseDDP<T>::nonlinear_rollout(T eps, HSDDP_OPTION &option)
{
    actual_cost = 0;
    max_pconstr = 0;
    max_tconstr = 0;
    DVec<T> xinit = x0; // initial condition for each phase
    DVec<T> xsim_init = x0;
    DVec<T> xend;
    DVec<T> xsim_end;
   
    for (int i = 0; i < n_phases; i++)
    {
        if (i > 0)
        {
            // If not the first phase, run resetmap at the end of previous phase
            phases[i - 1]->get_terminal_state(xend, xsim_end);
            xinit = phases[i - 1]->resetmap(xend);
            xsim_init = phases[i - 1]->resetmap(xsim_end);
        }

        phases[i]->set_initial_condition(xinit, xsim_init); // Set initial condition of current phase
        phases[i]->nonlinear_rollout(eps, option);             // run hybrid rollout for each phase
        
        // update the maximum constraint violations
        max_pconstr = std::min(max_pconstr, phases[i]->get_max_pconstrs()); // should have non-positive value
        max_tconstr = std::max(max_tconstr, phases[i]->get_max_tconstrs()); // should have non-negative value
    }
    
}

template <typename T>
bool MultiPhaseDDP<T>::line_search(HSDDP_OPTION &option)
{
    T exp_cost_change = 0;
    T exp_merit_change = 0;
    T eps = 1;
    T cost_prev = actual_cost;
    T merit_prev = merit;
    T feas_prev = feas;
    
    bool success = false;

#ifdef TIME_BENCHMARK
    fit_iter = 0;
#endif
    while (eps > 1e-5)
    {
        // nonlinear_rollout(eps, option);
        hybrid_rollout(eps, option);        
        compute_cost(option);
        feas = measure_dynamics_feasibility();
        merit = actual_cost + option.merit_rho * feas;

        exp_cost_change = eps * dV_1 + 0.5 * eps * eps * dV_2;
        exp_merit_change = exp_cost_change - eps * feas_prev;

#ifdef DEBUG
        printf("\t eps=%.3e \t actual change in cost=%.3e \t expeced change in cost=%.3e\n",
               eps, actual_cost - cost_prev, option.gamma * exp_cost_change);
#endif
        if (merit <= merit_prev + option.gamma * exp_merit_change)
        {
            success = true;
            break;
        }
        eps *= option.alpha;
#ifdef TIME_BENCHMARK
        fit_iter++;
#endif
    }
    return success;
}

template <typename T>
bool MultiPhaseDDP<T>::backward_sweep_regularized(T &regularization, HSDDP_OPTION &option)
{
    bool success = false;
    int iter = 1;

#ifdef TIME_BENCHMARK
    auto start = high_resolution_clock::now();
#endif

    while (!success)
    {
#ifdef DEBUG
        printf("\t regularization = %f\n", regularization);
#endif
        success = backward_sweep(regularization);
        if (success)
            break;

        regularization = std::max(regularization * option.update_regularization, T(1e-03));
        iter++;

        if (regularization > 1e2)
        {
            printf("Regularization term exceeds maximum value! \n");
            printf("Optimization terminates! \n");
            break;
        }
    }
#ifdef TIME_BENCHMARK
    stop = high_resolution_clock::now();
    duration = duration_ms(stop - start);
    time_per_iter.n_bws = iter;
    time_per_iter.time_bws = duration.count();
    time_per_iter.time_partial = time_partial;
    time_per_iter.DDP_iter = iter;
#endif
    regularization = regularization / 20;
    if (regularization < 1e-06)
        regularization = 0;
    return success;
}

/*
  @brief: perform backward sweep for the multi-phase problem
  @params:
          regularization: regularization parameter
  @return: success of backward sweep
*/
template <typename T>
bool MultiPhaseDDP<T>::backward_sweep(T regularization)
{
    bool success = true;
    T dVprime = 0;
    DVec<T> Gprime;
    DMat<T> Hprime;
    DVec<T> xend;
    DMat<T> Px;
    size_t xs(0), xs_next(0);
    // T dV_1_temp(0), dV_2_temp(0);
    // dV_1 = 0;
    // dV_2 = 0;
    for (int i = n_phases - 1; i >= 0; i--)
    {
        xs = phases[i]->get_state_dim();
        xs_next = (i < n_phases - 1) ? phases[i + 1]->get_state_dim() : xs;
        phases[i]->get_terminal_state(xend);
        Gprime.setZero(xs_next);
        Hprime.setZero(xs_next, xs_next);
        Px.setZero(xs, xs_next);
        // If not the last phase, udpate approximated value function back propagated from next phase
        if (i <= n_phases - 2)
        {
            phases[i]->resetmap_partial(Px, xend);
            phases[i + 1]->get_value_approx(Gprime, Hprime); // get the value infomation at the begining of phase i+1
            impact_aware_step(Gprime, Hprime, Px);
        }

        // If backward sweep of current phase fails, break and return false
        success = phases[i]->backward_sweep(regularization, Gprime, Hprime);
        if (!success)
            return success;

        // phases[i]->get_exp_cost_change(dV_1_temp, dV_2_temp);
        // dV_1 += dV_1_temp;
        // dV_2 += dV_2_temp;
    }
    return success;
}

template <typename T>
void MultiPhaseDDP<T>::solve(HSDDP_OPTION option)
{
    int iter = 0;
    int iter_ou = 0;
    int iter_in = 0;

    T cost_prev(0), merit_prev(0);
    bool success = false; // currently defined only for backward sweep
    bool ReB_active = option.ReB_active;    

#ifdef TIME_BENCHMARK
    time_ddp.clear();
    double time_partial = 0;
    TIME_PER_ITERATION time_per_iter;
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
#endif        
    
    hybrid_rollout(0, option); 
    update_nominal_trajectory();

    while (iter_ou < option.max_AL_iter)
    {
        iter_ou++;
        max_tconstr_prev = max_tconstr;
        max_pconstr_prev = max_pconstr;

#ifdef DEBUG
        printf("outer loop iteration %d \n", iter_ou);
#endif

#ifdef TIME_BENCHMARK
        start = high_resolution_clock::now();
#endif  
                  

#ifdef TIME_BENCHMARK
        stop = high_resolution_clock::now();
        duration = duration_ms(stop - start);
        time_partial = duration.count();
#endif
        
        T regularization = 0;
        iter_in = 0;
        while (iter_in < option.max_DDP_iter)
        {
            compute_cost(option);
            feas = measure_dynamics_feasibility();
            merit = actual_cost + option.merit_rho * feas;

            printf("total cost = %f, feasibility = %f \n", actual_cost, feas);
            iter_in++;
            iter++;
#ifdef DEBUG
            printf("\t inner loop iteration %d \n", iter_in);
#endif
            cost_prev = actual_cost;
            merit_prev = merit;

#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif
            LQ_approximation(option);
            success = backward_sweep_regularized(regularization, option);
            if (!success)
            {
                goto bad_solve;
            }
            
            if (option.MS)
            {
                linear_rollout(1.0, option);
            }
                        
            if (line_search(option))
            {
                // if line search succeeds, accept the step
                update_nominal_trajectory();                
            }
            else
            {
                // else do not update
                actual_cost = cost_prev;
                merit = merit_prev;
            }


#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_per_iter.n_fit = fit_iter;
            time_per_iter.time_fit = duration.count();
            time_ddp.push_back(time_per_iter);
#endif
           
            if ((fabs(cost_prev - actual_cost) < option.cost_thresh)&&(feas <= option.dynamics_feas_thresh))
                break;

#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif
            

#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_partial = duration.count();
#endif
        }
        if (option.AL_active)
        {
            update_AL_params(option);
        }
        if (option.ReB_active)
        {
            update_REB_params(option);
        }

#ifdef DEBUG
        printf("terminal constraint violation = %f \n", max_tconstr);
        printf("path constraint violation = %f \n", fabs(max_pconstr));
#endif

        if (max_tconstr < option.tconstr_thresh && fabs(max_pconstr) < option.pconstr_thresh && feas <= option.dynamics_feas_thresh)
        {
            printf("all constraints satisfied \n");
            printf("optimization finished \n");
            break;
        }
        if (fabs(max_tconstr - max_tconstr_prev) < 0.0001 && fabs(max_pconstr - max_pconstr_prev) < 0.0001 && feas <= option.dynamics_feas_thresh)
        {
            printf("constraints stop decreasing \n");
            printf("optimization terminates \n");
            break;
        }
    }
    if (iter_ou >= option.max_AL_iter)
    {
        printf("maximum iteration reached \n");
    }

    printf("total cost = %f \n", actual_cost);
    printf("terminal constraint violation = %f \n", max_tconstr);
    printf("path constraint violation = %f \n", fabs(max_pconstr));
    printf("dynamics feas = %f \n", feas);

bad_solve:
    if (!success)
    {
        printf(RED);
        printf("Failed to solve the optimization due to too large regularization \n");
        printf(RESET);
    }
}

template <typename T>
void MultiPhaseDDP<T>::compute_cost(const HSDDP_OPTION &option)
{
    actual_cost = 0;
    for (int i = 0; i < n_phases; i++)
    {
        phases[i]->compute_cost(option);
        actual_cost += phases[i]->get_actual_cost();
    }
}

template <typename T>
void MultiPhaseDDP<T>::LQ_approximation(HSDDP_OPTION &option)
{
    for (int i = 0; i < n_phases; i++)
    {
        phases[i]->LQ_approximation(option);
    }
}
/*
    @brief: Set a method to initialize the dynamics everytime berfore running forward sweep
    @Note:  In future version, this function need modified for general purpose. Currently it supports foothold
            initialization.
*/
template <typename T>
void MultiPhaseDDP<T>::set_dynamics_init_callback(function<void(DVec<T>)> dynamics_init_callback_)
{
    dynamics_init_callback = dynamics_init_callback_;
}

/*
    @brief: Always run before forward sweep to initialize the dyanmics model
*/
template <typename T>
void MultiPhaseDDP<T>::run_before_forward_sweep()
{
    if (nullptr != dynamics_init_callback)
    {
        dynamics_init_callback(x0);
    }
}

/*
  @brief: perform impact-aware backward DDP step
  @params:
          Px: reset map
          G, H: gradient and hessian of value function propaged from next phase through phase transition
*/

template <typename T>
void MultiPhaseDDP<T>::impact_aware_step(DVec<T> &G, DMat<T> &H, const DMat<T> &Px)
{
    G = Px.transpose() * G;
    H = Px.transpose() * H * Px;
}

template <typename T>
void MultiPhaseDDP<T>::update_REB_params(HSDDP_OPTION &option)
{
    for (auto &phase : phases)
    {
        phase->update_REB_params(option);
    }
}

template <typename T>
void MultiPhaseDDP<T>::update_AL_params(HSDDP_OPTION &option)
{
    for (auto &phase : phases)
    {
        phase->update_AL_params(option);
    }
}

template <typename T>
void MultiPhaseDDP<T>::update_nominal_trajectory()
{
    for (auto &phase : phases)
    {
        phase->update_nominal_trajectory();
    }
}

template <typename T>
T MultiPhaseDDP<T>::measure_dynamics_feasibility(int norm_id)
{
    T feasibility = 0;

    for (auto &phase : phases)
    {
        feasibility += phase->measure_dynamics_feasibility(norm_id);
    }

    if (norm_id == 2)
    {
        feasibility = sqrt(feasibility);
    }
    
    return feasibility;
}
template class MultiPhaseDDP<double>;
