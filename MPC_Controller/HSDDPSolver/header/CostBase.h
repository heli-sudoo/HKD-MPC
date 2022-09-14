#ifndef COSTBASE_H
#define COSTBASE_H

#include "HSDDP_CPPTypes.h"
#include <memory>
#include <vector>
#include <deque>

using std::vector;
using std::deque;
using std::shared_ptr;

// Running cost data structure
template<typename T, size_t xs_, size_t us_, size_t ys_>
struct RCostData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T l;
    VecM<T, xs_> lx;
    VecM<T, us_> lu;
    VecM<T, ys_> ly;
    MatMN<T, xs_, xs_> lxx;
    MatMN<T, us_, xs_> lux;
    MatMN<T, us_, us_> luu;
    MatMN<T, ys_, ys_> lyy;

    RCostData(){Zeros();}
    void Zeros()
    {
        l = 0;
        lx.setZero();
        lu.setZero();
        ly.setZero();
        lxx.setZero();
        luu.setZero();
        lux.setZero();
        lyy.setZero();
    }
    void add(const RCostData<T,xs_,us_,ys_>& cost_to_add)
    {
        l += cost_to_add.l;
        lx += cost_to_add.lx;
        lu += cost_to_add.lu;
        ly += cost_to_add.ly;
        luu += cost_to_add.luu;
        lux += cost_to_add.lux;
        lyy += cost_to_add.lyy;
    }
};

// Terminal cost data structure
template<typename T, size_t xs_>
struct TCostData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T Phi;
    VecM<T, xs_> Phix;
    MatMN<T, xs_, xs_> Phixx;

    TCostData(){Zeros();}
    void Zeros()
    {
        Phi = 0;
        Phix.setZero();
        Phixx.setZero();
    }
    void add(const TCostData<T,xs_>& cost_to_add)
    {
        Phi += cost_to_add.Phi;
        Phix += cost_to_add.Phix;
        Phixx += cost_to_add.Phixx;
    }
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
    virtual void running_cost(RCost&, const State& x, const  Contrl& u, const Output& y, T dt, int k) {}

    virtual void running_cost_par(RCost&, const State& x, const  Contrl& u, const Output& y, T dt, int k) {}

    virtual void terminal_cost(TCost&, const State& x, int kend) {}

    virtual void terminal_cost_par(TCost&, const State&x, int kend) {}

    virtual void penalize_foot_change(int l) {}

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

public:
    QuadraticCost();
    
    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void terminal_cost(TCost&, const State& x, int kend=0) override;

    void terminal_cost_par(TCost&, const State&x, int kend=0) override;

    void set_reference(deque<VecM<T, xs_>>* Xr_in,
                       deque<VecM<T, us_>>* Ur_in,
                       deque<VecM<T, ys_>>* Yr_in);

protected:
    MatMN<T, xs_, xs_> Q;
    MatMN<T, us_, us_> R;
    MatMN<T, ys_, ys_> S;
    MatMN<T, xs_, xs_> Qf;

private:
     // iterators to state, control, output trajectories
    deque<VecM<T, xs_>>* Xr = nullptr;
    deque<VecM<T, us_>>* Ur = nullptr;
    deque<VecM<T, ys_>>* Yr = nullptr;
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
    vector<shared_ptr<CostBase<T, xs_, us_, ys_>>> cost_ptrs;

public:
    void add_cost(shared_ptr<CostBase<T, xs_, us_, ys_>> ptr_to_cost_to_add){
        cost_ptrs.push_back(ptr_to_cost_to_add);
    }
    
    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0);

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0);

    void terminal_cost(TCost&, const State& x, int kend=0);

    void terminal_cost_par(TCost&, const State&x, int kend=0);

private:
    RCostData<T, xs_, us_, ys_> rcostdata_temp;
    TCostData<T, xs_> tcostdata_temp;
};  
#endif // COSTBASE_H