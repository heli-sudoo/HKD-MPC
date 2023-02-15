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

    virtual void set_initial_condition(DVec<T>& x0_, DVec<T>& xsim_0_);

    virtual void set_initial_condition_dx(DVec<T>& dx0_) = 0; // initial condition for linear roll-out

    virtual void set_nominal_initial_condition(DVec<T>& x0_){} // initial condition of the nominal trajectory

    virtual void forward_sweep(T eps, HSDDP_OPTION&, bool) = 0; // compute both dynamics propagation and dynamics partials

    virtual void linear_rollout(T eps, HSDDP_OPTION&) = 0;

    virtual void hybrid_rollout(T eps, HSDDP_OPTION&) = 0;

    virtual void nonlinear_rollout(T eps, HSDDP_OPTION&) = 0;

    virtual void LQ_approximation(HSDDP_OPTION&) = 0;

    virtual bool backward_sweep(T regularization, T dVprime, DVec<T> Gprime, DMat<T> Hprime) = 0;
    
    virtual DVec<T> resetmap(DVec<T>&) = 0;

    virtual void resetmap_partial(DMat<T>&, DVec<T>&) = 0;

    virtual void get_value_info_at_init(T& dV, DVec<T>& G, DMat<T>& H) = 0;

    virtual DVec<T> get_terminal_state() = 0;

    virtual DVec<T> get_terminal_state_dx() = 0;

    virtual DVec<T> get_terminal_state_nominal(){}

    virtual T get_actual_cost() = 0;   

    virtual T get_max_tconstrs() {}

    virtual T get_max_pconstrs() {}

    virtual size_t get_state_dim() {return 0;}

    virtual size_t get_control_dim() {return 0;}

    virtual void update_AL_params(HSDDP_OPTION& ) {} // do nothing

    virtual void update_REB_params(HSDDP_OPTION& ) {} // do nothing

    virtual void update_nominal_trajectory() = 0;

    virtual void empty_control(){}

    virtual void push_back() {}

    virtual void pop_front() {}

    virtual void reset_params() {}

    virtual T check_dynamics_feasibility() {return 0;}

};

#endif //SINGLEPHASE_BASE_H