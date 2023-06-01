#ifndef SINGLEPHASE_HSDDP
#define SINGLEPHASE_HSDDP

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <memory>
#include <functional>
#include <numeric>
#include "HSDDP_CompoundTypes.h"
#include "SinglePhaseBase.h"
#include "TrajectoryManagement.h"
#include "ConstraintsBase.h"
#include "SinglePhaseInterface.h"

using std::string;
using std::deque;
using std::shared_ptr;
using std::function;

template<typename> class MultiPhaseDDP; //forward declaration of MultiPhaseDDP class template

template<typename T, size_t xs, size_t us, size_t ys>
class SinglePhase: public SinglePhaseBase<T>
{
public:
    typedef VecM<T, xs> State;
    typedef VecM<T, us> Contrl;
    typedef VecM<T, ys> Output;
    typedef MatMN<T, xs, xs> StateMap;
    typedef MatMN<T, xs, us> ContrlMap;
    typedef MatMN<T, ys, xs> OutputMap;
    typedef MatMN<T, ys, us> DirectMap;

    friend class MultiPhaseDDP<T>;

private:
    // function wrapper of Callable dynamics
    function<void(State&, Output&, State&, Contrl&, T)> dynamics;

    // function wrapper of Callable dynamics linearizaiton
    function<void(StateMap&, ContrlMap&, OutputMap&, DirectMap&, State&, Contrl&, T)> dynamics_partial;

    // function wrapper of Callable resetmap
    function<void(DVec<T>&, DVec<T>&)> resetmap_func_handle;

    // function wrapper of Callable resetmap partial
    function<void(DMat<T>&, DVec<T>&)> resetmap_partial_func_handle;

    // constraint Contaniner
    ConstraintContainer<T, xs, us, ys> constraintContainer;

    // pointer to cost object
    CostContainer<T,xs,us,ys> costContainer;


public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SinglePhase();
    
    void set_trajectory(shared_ptr<Trajectory<T,xs,us,ys>> traj_);

    void set_dynamics(function<void(State&, Output&, State&, Contrl&, T)> dynamics_in){
        dynamics = dynamics_in;
    }

    void set_dynamics_partial(function<void(StateMap&, ContrlMap&, OutputMap&, DirectMap&, State&, Contrl&, T)> dynamics_partial_in){
        dynamics_partial = dynamics_partial_in;
    }

    void set_resetmap(function<void(DVec<T>&, DVec<T>&)> resetmap_func_handle_in){
        resetmap_func_handle = resetmap_func_handle_in;
    }

    void set_resetmap_partial(function<void(DMat<T>&, DVec<T>&)> resetmap_partial_in){
        resetmap_partial_func_handle = resetmap_partial_in;
    }

    void set_time_offset(float t_offset_in) {t_offset = t_offset_in;}

    void add_cost(shared_ptr<CostBase<T,xs,us,ys>> ptr_to_cost_to_add){
        costContainer.add_cost(ptr_to_cost_to_add);
    }
    
    void add_pathConstraint(shared_ptr<PathConstraintBase<T, xs, us, ys>> pConstraint){
        constraintContainer.add_pathConstraint(pConstraint);
    }

    void add_terminalConstraint(shared_ptr<TerminalConstraintBase<T, xs>> tConstraint){
        constraintContainer.add_terminalConstraint(tConstraint);
    }    

public:
    void initialization() override;

    void set_initial_condition(DVec<T>& x0_) override;

    void set_initial_condition(DVec<T>& x0_, DVec<T>& x_sim_0_) override; // set initial condition for both local states and simulation states

    void set_initial_condition_dx(DVec<T>& dx0_) override; // set initial condition for linear rollout

    void set_nominal_initial_condition(DVec<T>& x0_) override;

    void warmstart() override {}    

    void linear_rollout(T eps, HSDDP_OPTION&) override;

    bool hybrid_rollout(T eps, HSDDP_OPTION&, bool is_last_phase = false) override;

    void LQ_approximation(HSDDP_OPTION&) override;

    bool backward_sweep(T regularization, DVec<T> Gprime, DMat<T> Hprime) override;

    DVec<T> resetmap(DVec<T>& x_) override;

    void resetmap_partial(DMat<T>& Px_, DVec<T>& x_) override;    

    void get_value_approx(DVec<T>& G_out, DMat<T>& H_out) override;

    void get_exp_cost_change(T& dV_1_out, T& dV_2_out) override;

    void get_terminal_state(DVec<T>& xend_out) override;

    void get_terminal_state(DVec<T>& xend_out, DVec<T>& xsim_end_out) override;

    void get_terminal_state_dx(DVec<T>& dx_end_out) override;

    T get_actual_cost() override;

    T get_max_tconstrs() override;

    T get_max_pconstrs() override;

    size_t get_state_dim() override {return xs;}

    size_t get_control_dim() override {return us;}

    void update_AL_params(HSDDP_OPTION& option) override;

    void update_REB_params(HSDDP_OPTION& option) override;

    void update_nominal_trajectory() override;

    void empty_control() override{
        traj->zero_val_approx();
    }

    void push_back_default() override;

    void pop_front() override;

    void reset_params() override {constraintContainer.reset_params();}

    T measure_dynamics_feasibility(int norm_id=2) override {return traj->measure_dynamics_feasibility(norm_id);}

    void compute_defect() {traj->compute_defect();}

    void update_SS_config(int ss_sz) override {
        SS_set.resize(ss_sz);        
        std::iota(SS_set.begin(), SS_set.end(), 0);
        }

    void get_trajectory(std::vector< std::vector<float>> & x_tau, 
                        std::vector< std::vector<float>> & u_tau) override ;    

    void print() override;                         

private:
    void update_trajectory_ptrs();    

    void compute_cost(const HSDDP_OPTION& option);                                                                                                

    void update_running_cost_with_pconstr(size_t k);

    void update_running_cost_par_with_pconstr(size_t k);
                                                                                        
    void update_running_cost_with_smooth();    

    void update_terminal_cost_with_tconstr();           

    void update_terminal_cost_par_with_tconstr(); 
                                                                                                                         

private:
    shared_ptr<Trajectory<T,xs,us,ys>> traj;    // pointer to single phase trajectory
    int phase_horizon;                          // number of timesteps for one phase (number of states = phase_horizon + 1)
    T dt;                                       // simulation timestep
    float t_offset;                             // time offset w.r.t. to the begining of the entire multi-phase trajectory

   /* pointers to hold state, control and output trajectory */    
    deque<VecM<T, xs>>* Xbar = nullptr;
    deque<VecM<T, xs>>* X = nullptr;
    deque<VecM<T, us>>* Ubar = nullptr;
    deque<VecM<T, us>>* U = nullptr;
    deque<VecM<T, ys>>* Y = nullptr;

    /* pointers to hold simulation state, defect etc */
    deque<VecM<T, xs>> *Xsim = nullptr;
    deque<VecM<T, xs>> *Defect = nullptr;
    deque<VecM<T, xs>> *Defect_bar = nullptr;

    /* pointers to hold linearized dynamics */
    deque<MatMN<T, xs, xs>>* A = nullptr;
    deque<MatMN<T, xs, us>>* B = nullptr;
    deque<MatMN<T, ys, xs>>* C = nullptr;
    deque<MatMN<T, ys, us>>* D = nullptr;

    /* pointers to hold value function approximation */
    deque<T>* V = nullptr;
    deque<T>* dV = nullptr;
    deque<VecM<T, us>>* dU = nullptr;
    deque<VecM<T, xs>>* G = nullptr;
    deque<MatMN<T, xs, xs>>* H = nullptr;
    deque<MatMN<T, us, xs>>* K = nullptr;
    
    deque<VecM<T, xs>>* dX = nullptr;
    

    // running cost and terminal cost
    deque<RCostData<T, xs, us, ys>>* rcostData = nullptr;
    TCostData<T, xs>* tcostData  = nullptr;             
    
    // Q matrices 
    VecM<T, xs> Qx;
    VecM<T, us> Qu;
    MatMN<T, xs, xs> Qxx;
    MatMN<T, us, us> Quu;
    MatMN<T, us, xs> Qux;
    MatMN<T, xs, xs> Ixx;
    MatMN<T, us, us> Iuu;   
    Chol<T> Qxx_chol;
    Chol<T> Quu_chol;
    MatMN<T, us, us> Quu_inv;

    size_t tConstrs_size = 0;
    size_t pConstrs_size = 0;
    DVec<T> Bar;
    DVec<T> Bard;
    DVec<T> Bardd;

    // Initial condition
    VecM<T, xs> x_init;
    VecM<T, xs> xsim_init;
    VecM<T, xs> dx_init;    

    // Shooting state set
    std::vector<int> SS_set;

    T actual_cost;
    T dV_1;     // expected cost change (first-order term)
    T dV_2;     // expected cost change (second-order term)
};

#endif // SINGLEPHASE_HSDDP