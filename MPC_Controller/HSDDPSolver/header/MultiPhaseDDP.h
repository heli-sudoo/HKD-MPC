#ifndef MULTIPHASEDDP_H
#define MULTIPHASEDDP_H

#include <vector>
#include <deque>
#include <memory> // smart pointer (since C++ 11)
#include <functional>
#include "SinglePhase.h"

using std::vector;
using std::deque;
using std::shared_ptr;
using std::function;

template <typename T>
class MultiPhaseDDP
{
private:
    deque<shared_ptr<SinglePhaseBase<T>>> phases;

public:
    MultiPhaseDDP() {} // default constructor

    void set_multiPhaseProblem(deque<shared_ptr<SinglePhaseBase<T>>> phases_in){
        phases = phases_in;
        n_phases = phases.size();
        actual_cost = 0;        
        max_pconstr = 0;
        max_pconstr_prev = 0;
        max_tconstr = 0;
        max_tconstr_prev = 0;
    }

    void set_initial_condition(DVec<T> x0_in) { x0 = x0_in; dx0.setZero(x0.size()); }

    void solve(HSDDP_OPTION option);

    void set_dynamics_init_callback(function<void(DVec<T>)> dynamics_init_callback_);

public:
    void forward_sweep(T eps, HSDDP_OPTION& option, int);

    void linear_rollout(T eps, HSDDP_OPTION& option);

    void hybrid_rollout(T eps, HSDDP_OPTION& option);

    bool line_search(HSDDP_OPTION& option);

    void LQ_approximation(HSDDP_OPTION& option);

    bool backward_sweep(T regularization);

    bool backward_sweep_regularized(T& regularization, HSDDP_OPTION&option);

    bool forward_iteration(HSDDP_OPTION& option);   

    void impact_aware_step(DVec<T>&G, DMat<T>&H, const DMat<T>& Px);

    void update_nominal_trajectory();

    void run_before_forward_sweep();

    void update_AL_params(HSDDP_OPTION& option);

    void update_REB_params(HSDDP_OPTION& option);

    T get_actual_cost() {return actual_cost;}

    // T get_exp_cost_change() {return exp_cost_change;}

    void empty_solution(){
        for (auto phase:phases)
        {
            phase->empty_control();
        }        
        phases.clear();
        actual_cost = 0;
        // exp_cost_change = 0;
        max_tconstr = 0;
        max_pconstr = 0;
        max_tconstr_prev = 0;
        max_pconstr_prev = 0;
    }

    T check_dynamics_feasibility();    

private:
    int n_phases;

    T actual_cost;
    // T exp_cost_change;
    T dV_1;         // expected cost change due to first order
    T dV_2;         // expected cost change due to second order

    T max_tconstr_prev;
    T max_pconstr_prev;
    T max_tconstr; // maximum terminal contraint violation of all time
    T max_pconstr; // maximum path constraint violation of all time

    DVec<T> x0;
    DVec<T> dx0;
private:
    function<void(DVec<T>)> dynamics_init_callback;

};

#endif // MULTIPHASEDDP_H