#ifndef SINGLEPHASE_BASE_H
#define SINGLEPHASE_BASE_H

#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"


template<typename> class MultiPhaseDDP; //forward declaration of MultiPhaseDDP class template

template<typename T>
class SinglePhaseBase
{

private:    
    friend class MultiPhaseDDP<T>; 

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SinglePhaseBase(){}
                            
    virtual void warmstart() = 0;

    virtual void initialization() = 0;

    virtual void set_initial_condition(DVec<T>& x0_) = 0; // initial condition of the acutual trajectory

    virtual void set_initial_condition(DVec<T>& x0_, DVec<T>& xsim_0_) {(void) (x0_); (void) (xsim_0_);}

    virtual void set_initial_condition_dx(DVec<T>& dx0_) = 0; // initial condition for linear roll-out

    virtual void set_nominal_initial_condition(DVec<T>& x0_){(void) (x0_);} // initial condition of the nominal trajectory

    virtual void linear_rollout(T eps, HSDDP_OPTION&) = 0;

    virtual bool hybrid_rollout(T eps, HSDDP_OPTION&, bool is_last_phase = false) = 0;    

    virtual void LQ_approximation(HSDDP_OPTION&) = 0;

    virtual bool backward_sweep(T regularization, DVec<T> Gprime, DMat<T> Hprime) = 0;
    
    virtual DVec<T> resetmap(DVec<T>&) = 0;

    virtual void resetmap_partial(DMat<T>& Px, DVec<T>& x) = 0;

    virtual void get_value_approx(DVec<T>& G, DMat<T>& H) = 0;

    virtual void get_exp_cost_change(T& dV_1, T& dV_2) = 0;

    virtual void get_terminal_state(DVec<T>& xend) = 0;

    virtual void get_terminal_state(DVec<T>& xend, DVec<T>& xsim_end) = 0;

    virtual void get_terminal_state_dx(DVec<T>& dx_end) = 0;

    virtual T get_actual_cost() = 0;   

    virtual T get_max_tconstrs() {return (T)(0);}

    virtual T get_max_pconstrs() {return (T)(0);}

    virtual size_t get_state_dim() {return 0;}

    virtual size_t get_control_dim() {return 0;}

    virtual void update_AL_params(HSDDP_OPTION& ) {} // do nothing

    virtual void update_REB_params(HSDDP_OPTION& ) {} // do nothing

    virtual void update_nominal_trajectory() = 0;

    virtual void empty_control(){}

    virtual void push_back_default() {}

    virtual void pop_front() {}

    virtual void reset_params() {}

    virtual T measure_dynamics_feasibility(int norm_id) {(void) (norm_id); return 0;}    

    virtual void update_SS_config(int ss_sz) {(void) (ss_sz);}

    virtual void compute_cost(const HSDDP_OPTION& option) = 0;

    virtual void get_trajectory(std::vector< std::vector<float>> & x_tau, std::vector< std::vector<float>> & u_tau) 
    {(void) x_tau; (void) (u_tau);}

    virtual void print() {}

};

#endif //SINGLEPHASE_BASE_H