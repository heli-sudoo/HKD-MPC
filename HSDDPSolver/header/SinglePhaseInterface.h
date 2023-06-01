#ifndef COSTBASE_H
#define COSTBASE_H

#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"
#include <memory>       // smart pointer
#include <vector>
#include <string>

using std::shared_ptr;


template <size_t xs, size_t us, size_t ys>
class SinglePhaseReferenceAbstract
{
public:
    SinglePhaseReferenceAbstract(){}

    virtual void get_reference_at_t (VecM<double, xs>& xr, VecM<double, us>& ur, VecM<double, ys>& yr, float t) {
        (void) (xr);
        (void) (ur);
        (void) (yr);
        (void) (t);
        printf("Need to wrap around the function get_reference_at_t for your problem \n");
    };

    virtual void get_reference_at_t (VecM<double, xs>& xr, float t){
        (void) (xr);
        (void) (t);
        printf("Need to wrap around the function get_reference_at_t for your problem \n");
    };
};

template <typename T, size_t xs_, size_t us_, size_t ys_>
class CostBase 
{
public:    
    typedef VecM<T, xs_> State;
    typedef VecM<T, us_> Contrl;
    typedef VecM<T, ys_> Output;
    typedef RCostData<T,xs_,us_,ys_> RCost;
    typedef TCostData<T,xs_> TCost;

public:    
    CostBase(const std::string& cost_name_in) {cost_name = cost_name_in;}

    virtual void running_cost(RCost&, const State& x, const  Contrl& u, const Output& y, T dt, float t) = 0;

    virtual void running_cost_par(RCost&, const State& x, const  Contrl& u, const Output& y, T dt, float t) = 0;

    virtual void terminal_cost(TCost&, const State& x, float tend) = 0;

    virtual void terminal_cost_par(TCost&, const State&x, float tend) = 0;

public:
    std::string cost_name;    
};

template <typename T, size_t xs_, size_t us_, size_t ys_>
class QuadraticCost:public CostBase<T,xs_,us_,ys_>
{
public:
    using typename CostBase<T,xs_,us_,ys_>::State;
    using typename CostBase<T,xs_,us_,ys_>::Contrl;
    using typename CostBase<T,xs_,us_,ys_>::Output;
    using typename CostBase<T,xs_,us_,ys_>::RCost;
    using typename CostBase<T,xs_,us_,ys_>::TCost;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    QuadraticCost(const std::string& cost_name_in):
        CostBase<T,xs_,us_,ys_>(cost_name_in){
        /* default weightings */
        Q.setIdentity();
        R.setIdentity();
        S.setZero();
        Qf.setIdentity();
    }

    QuadraticCost() : CostBase<T,xs_,us_,ys_>("Quadratic Cost"){
        /* default weightings */
        Q.setIdentity();
        R.setIdentity();
        S.setZero();
        Qf.setIdentity();
    }

    virtual void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0) override;

    virtual void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0) override;

    virtual void terminal_cost(TCost&, const State& x, float tend = 0) override;

    virtual void terminal_cost_par(TCost&, const State&x, float tend = 0) override;

protected:
    MatMN<T, xs_, xs_> Q;
    MatMN<T, us_, us_> R;
    MatMN<T, ys_, ys_> S;
    MatMN<T, xs_, xs_> Qf;
};

template <typename T, size_t xs_, size_t us_, size_t ys_>
class QuadraticTrackingCost: public QuadraticCost<T, xs_, us_, ys_>
{
public:
    using typename CostBase<T,xs_,us_,ys_>::State;
    using typename CostBase<T,xs_,us_,ys_>::Contrl;
    using typename CostBase<T,xs_,us_,ys_>::Output;
    using typename CostBase<T,xs_,us_,ys_>::RCost;
    using typename CostBase<T,xs_,us_,ys_>::TCost;

public: 
    QuadraticTrackingCost(const std::string& cost_name_in):
    QuadraticCost<T, xs_, us_, ys_>(cost_name_in)
    {
        xr_t.setZero();
        ur_t.setZero();
        yr_t.setZero();
    }    

    QuadraticTrackingCost()
    {
        xr_t.setZero();
        ur_t.setZero();
        yr_t.setZero();
    }

    virtual void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0) override;

    virtual void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0) override;

    virtual void terminal_cost(TCost&, const State& x, float tend = 0) override;

    virtual void terminal_cost_par(TCost&, const State&x, float tend = 0) override;

    void set_reference(SinglePhaseReferenceAbstract<xs_, us_, ys_>* reference_in) {reference = reference_in;}

private:
    SinglePhaseReferenceAbstract<xs_, us_, ys_>* reference = nullptr;

    VecM<double, xs_> xr_t;
    VecM<double, us_> ur_t;
    VecM<double, ys_> yr_t;
};

template <typename T, size_t xs_, size_t us_, size_t ys_>
class CostContainer
{
public:
    typedef VecM<T, xs_> State;
    typedef VecM<T, us_> Contrl;
    typedef VecM<T, ys_> Output;
    typedef RCostData<T,xs_,us_,ys_> RCost;
    typedef TCostData<T,xs_> TCost;

public:
    std::vector<shared_ptr<CostBase<T, xs_, us_, ys_>>> cost_ptrs;

public:
    void add_cost(shared_ptr<CostBase<T, xs_, us_, ys_>> ptr_to_cost_to_add){
        cost_ptrs.push_back(ptr_to_cost_to_add);
    }
    
    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0);

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t = 0);

    void terminal_cost(TCost&, const State& x, float tend = 0);

    void terminal_cost_par(TCost&, const State&x, float tend = 0);

private:
    RCostData<T, xs_, us_, ys_> rcostdata_temp;
    TCostData<T, xs_> tcostdata_temp;
};  
#endif // COSTBASE_H