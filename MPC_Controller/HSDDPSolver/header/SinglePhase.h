#ifndef SINGLEPHASE_HSDDP
#define SINGLEPHASE_HSDDP

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <memory>
#include <functional>
#include "HSDDP_CompoundTypes.h"
#include "SinglePhaseBase.h"
#include "TrajectoryManagement.h"
#include "ConstraintsBase.h"
#include "CostBase.h"

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

public:
    friend class MultiPhaseDDP<T>;
    // function wrapper of Callable dynamics
    function<void(State&, Output&, State&, Contrl&)> dynamics;
    // function wrapper of Callable dynamics linearizaiton
    function<void(StateMap&, ContrlMap&, OutputMap&, DirectMap&, State&, Contrl&)> dynamics_partial;
    // function wrapper of Callable resetmap
    function<void(DVec<T>&, DVec<T>&)> resetmap_func_handle;
    // function wrapper of Callable resetmap partial
    function<void(DMat<T>&, DVec<T>&)> resetmap_partial_func_handle;

    // constraint Contaniner
    ConstraintContainer<T, xs, us, ys> constraintContainer;
    // pointer to cost object
    // shared_ptr<CostBase<T,xs, us, ys>> costContainer;    
    CostContainer<T,xs,us,ys> costContainer;


public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SinglePhase();
    
    void set_trajectory(shared_ptr<Trajectory<T,xs,us,ys>> traj_);

    void set_dynamics(function<void(State&, Output&, State&, Contrl&)> dynamics_in){
        dynamics = dynamics_in;
    }

    void set_dynamics_partial(function<void(StateMap&, ContrlMap&, OutputMap&, DirectMap&, State&, Contrl&)> dynamics_partial_in){
        dynamics_partial = dynamics_partial_in;
    }

    void set_resetmap(function<void(DVec<T>&, DVec<T>&)> resetmap_func_handle_in){
        resetmap_func_handle = resetmap_func_handle_in;
    }

    void set_resetmap_partial(function<void(DMat<T>&, DVec<T>&)> resetmap_partial_in){
        resetmap_partial_func_handle = resetmap_partial_in;
    }

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

    void forward_sweep(T eps, HSDDP_OPTION& option, bool) override;

    void linear_rollout(T eps, HSDDP_OPTION&) override;

    void hybrid_rollout(T eps, HSDDP_OPTION&) override;

    void nonlinear_rollout(T eps, HSDDP_OPTION&) override;

    void LQ_approximation(HSDDP_OPTION&) override;

    bool backward_sweep(T regularization, T dVprime, DVec<T> Gprime, DMat<T> Hprime) override;

    DVec<T> resetmap(DVec<T>& x_) override;

    void resetmap_partial(DMat<T>& Px_, DVec<T>& x_) override;

    void get_value_info_at_init(T&dV, DVec<T> &G, DMat<T> &H) override;

    DVec<T> get_terminal_state() override;

    DVec<T> get_terminal_state_dx() override;

    DVec<T> get_terminal_state_nominal() override;

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

    void push_back() override;

    void pop_front() override;

    void reset_params() override {constraintContainer.reset_params();}

    T check_dynamics_feasibility() override {return traj->dynamics_feasibility();}

private:
    void update_trajectory_ptrs();
    
    void compute_barrier(vector<IneqConstrData<T,xs,us,ys>>&, vector<REB_Param_Struct<T>>&); 
                                                                                            
    void update_running_cost_with_pconstr(RCostData<T,xs,us,ys>& rcost,
                                          vector<IneqConstrData<T,xs,us,ys>>& pconstrsData,
                                          vector<REB_Param_Struct<T>>& reb_params,
                                          bool flag);

    void update_running_cost_with_smooth();

    void update_terminal_cost_with_tconstr(vector<TConstrData<T, xs>>& tconstrsData,
                                           vector<AL_Param_Struct<T>>& al_params,
                                           bool flag);

private:
    shared_ptr<Trajectory<T,xs,us,ys>> traj;
    // BriefTrajectoryIterator<T,xs,us,ys> ref; // reference trajectory iterator
    int offset; // offset of the current phase to the start of the trajectory
    int phase_horizon;
    T dt;
    
   /* pointers to hold state, control and output trajectory */    
    deque<VecM<T, xs>>* Xbar = nullptr;
    deque<VecM<T, xs>>* X = nullptr;
    deque<VecM<T, us>>* Ubar = nullptr;
    deque<VecM<T, us>>* U = nullptr;
    deque<VecM<T, ys>>* Y = nullptr;

    /* pointers to hold simulation state, defect etc */
    deque<VecM<T, xs>> *Xsim = nullptr;
    deque<VecM<T, xs>> *Defect = nullptr;

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
    deque<MatMN<T, xs, xs>>* K = nullptr;
    
    deque<VecM<T, xs>>* dX = nullptr;
    

    // running cost and terminal cost
    deque<RCostData<T, xs, us, ys>>* rcostData = nullptr;
    TCostData<T, xs>* tcostData  = nullptr;  

    // constraints data


    // iterators to state, control, output trajectories
    typename deque<VecM<T, xs>>::iterator Xr;
    typename deque<VecM<T, us>>::iterator Ur;
    typename deque<VecM<T, ys>>::iterator Yr;
    
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

    T actual_cost;
};

#endif // SINGLEPHASE_HSDDP