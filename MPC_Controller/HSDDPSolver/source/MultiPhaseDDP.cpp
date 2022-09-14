#include "MultiPhaseDDP.h"
#include <algorithm>
#include "HSDDP_CompoundTypes.h"
#include "cTypes.h"

// #define DEBUG

#ifdef TIME_BENCHMARK
#include <chrono>
using namespace std::chrono;
using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;
static int fit_iter = 0;
#endif // TIME_BENCHMARK

/*
  @brief: perform forward sweep for the multi-phase problem
  @params:
          eps: line search parameter (0, 1)
*/
template <typename T>
void MultiPhaseDDP<T>::forward_sweep(T eps, HSDDP_OPTION &option, bool calc_partial)
{
    actual_cost = 0;
    max_pconstr = 0;
    max_tconstr = 0;
    DVec<T> xinit = x0; // initial condition for each phase

    run_before_forward_sweep();
    for (size_t i = 0; i < n_phases; i++)
    {
        if (i == 0)
        {
            phases[i]->set_nominal_initial_condition(x0);
        }

        if (i > 0)
        {
            // If not the first phase, run resetmap at the end of previous phase
            DVec<T> xend = phases[i - 1]->get_terminal_state();
            xinit = phases[i - 1]->resetmap(xend);
        }

        phases[i]->set_initial_condition(xinit);             // Set initial condition of current phase
        phases[i]->forward_sweep(eps, option, calc_partial); // run forward sweep for current phase
        actual_cost += phases[i]->get_actual_cost();         // update total cost
        // update the maximum constraint violations
        max_pconstr = std::min(max_pconstr, phases[i]->get_max_pconstrs()); // should have non-positive value
        max_tconstr = std::max(max_tconstr, phases[i]->get_max_tconstrs()); // should have non-negative value
    }
    // DVec<T> xend = phases.back()->get_terminal_state();
    // std::cout << "xend = " << xend.transpose() << std::endl;
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
    exp_cost_change = 0;
    T dVprime = 0;
    DVec<T> Gprime;
    DMat<T> Hprime;
    DVec<T> xend;
    DMat<T> Px;
    size_t xs(0), xs_next(0);
    for (int i = n_phases - 1; i >= 0; i--)
    {
        xs = phases[i]->get_state_dim();
        xs_next = (i < n_phases - 1) ? phases[i + 1]->get_state_dim() : xs;
        xend = phases[i]->get_terminal_state();
        Gprime.setZero(xs_next);
        Hprime.setZero(xs_next, xs_next);
        Px.setZero(xs, xs_next);
        // If not the last phase, udpate approximated value function back propagated from next phase
        if (i <= n_phases - 2)
        {
            phases[i]->resetmap_partial(Px, xend);
            phases[i + 1]->get_value_info_at_init(dVprime, Gprime, Hprime); // get the value infomation at the begining of phase i+1
            impact_aware_step(Gprime, Hprime, Px);
        }
        // If backward sweep of current phase fails, break and return false
        success = phases[i]->backward_sweep(regularization, dVprime, Gprime, Hprime);
        if (!success)
            return success;
    }

    T dV_0 = 0;
    DVec<T> G_0;
    DMat<T> H_0;
    phases[0]->get_value_info_at_init(dV_0, G_0, H_0);

    exp_cost_change = dV_0;
    return success;
}

template <typename T>
bool MultiPhaseDDP<T>::forward_iteration(HSDDP_OPTION &option)
{
    T eps = 1;
    T cost_prev = actual_cost;
    bool success = false;
#ifdef TIME_BENCHMARK
    fit_iter = 0;
#endif
    while (eps > 1e-5)
    {
        forward_sweep(eps, option, false); // perform forward sweep but not computing dynamics linearization
#ifdef DEBUG
        printf("\t eps=%.3e \t actual change in cost=%.3e \t expeced change in cost=%.3e\n",
               eps, actual_cost - cost_prev, option.gamma * eps * (1 - eps / 2) * exp_cost_change);
#endif
        if (actual_cost <= cost_prev + option.gamma * eps * (1 - eps / 2) * exp_cost_change)
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
void MultiPhaseDDP<T>::solve(HSDDP_OPTION option)
{
    int iter = 0;
    int iter_ou = 0;
    int iter_in = 0;

    T cost_prev;
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

    while (iter_ou < option.max_AL_iter)
    {
        iter_ou++;
        max_tconstr_prev = max_tconstr;
        max_pconstr_prev = max_pconstr;

#ifdef DEBUG
        printf("outer loop iteration %d \n", iter_ou);
#endif

        option.ReB_active = ReB_active;
        // if (max_tconstr > 0.1 || 1 == iter_ou)
        if (max_tconstr > 0.1)
        {
            /* If terminal constraint violation is too large, turn path constraint off */
            option.ReB_active = 0;
        }

#ifdef TIME_BENCHMARK
        start = high_resolution_clock::now();
#endif
        forward_sweep(0, option, true);    
        printf("total cost = %f \n", actual_cost);   

#ifdef TIME_BENCHMARK
        stop = high_resolution_clock::now();
        duration = duration_ms(stop - start);
        time_partial = duration.count();
#endif

        update_nominal_trajectory();
        T regularization = 0;
        iter_in = 0;
        while (iter_in < option.max_DDP_iter)
        {
            iter_in++;
            iter++;
#ifdef DEBUG
            printf("\t inner loop iteration %d \n", iter_in);
#endif
            cost_prev = actual_cost;

#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif
            success = backward_sweep_regularized(regularization, option);
            if (!success)
            {
                goto bad_solve;
            }
            
            if (forward_iteration(option))
            {
                // if line search succeeds
                // accept the step
                update_nominal_trajectory();
            }
            else
            {
                // else do not update
                actual_cost = cost_prev;
            }

#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_per_iter.n_fit = fit_iter;
            time_per_iter.time_fit = duration.count();
            time_ddp.push_back(time_per_iter);
#endif
            // If cost change small, accept the DDP solution
            if (cost_prev - actual_cost < option.DDP_thresh)
                break;

#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif

            forward_sweep(0, option, true);

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

        if (max_tconstr < option.tconstr_thresh && fabs(max_pconstr) < option.pconstr_thresh)
        {
            printf("all constraints satisfied \n");
            printf("optimization finished \n");
            break;
        }
        if (fabs(max_tconstr - max_tconstr_prev) < 0.0001 && fabs(max_pconstr - max_pconstr_prev) < 0.0001)
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

    bad_solve:
    if (!success)
    {
        printf(RED);
        printf("Failed to solve the optimization due to too large regularization \n");
        printf(RESET);
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

template class MultiPhaseDDP<double>;
