#include "CostBase.h"
#include "HSDDP_Utils.h"
template <typename T, size_t xs_, size_t us_, size_t ys_>
QuadraticCost<T, xs_, us_, ys_>::QuadraticCost()
{
    /* default weightings */
    Q.setIdentity();
    R.setIdentity();
    S.setZero();
    Qf.setIdentity();
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::set_reference(deque<VecM<T, xs_>>* Xr_in,
                                                    deque<VecM<T, us_>>* Ur_in,
                                                    deque<VecM<T, ys_>>* Yr_in)
{
    Xr = Xr_in;
    Ur = Ur_in;
    Yr = Yr_in;
} 

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::running_cost(RCost& rcost, const State& x, const Contrl& u, 
                                                   const Output& y, T dt, int k)
{    
    const auto& dx = x - Xr->at(k);
    const auto& du = u - Ur->at(k);
    const auto& dy = y - Yr->at(k);
    rcost.l = 0.5*dx.transpose() * Q * dx;
    rcost.l += 0.5*du.transpose() * R * du;
    rcost.l += 0.5*dy.transpose() * S * dy;
    rcost.l *= dt;
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::running_cost_par(RCost& rcost, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, int k)
{
    const auto& dx = x - Xr->at(k);
    const auto& du = u - Ur->at(k);
    const auto& dy = y - Yr->at(k);    
    rcost.lx = dt * Q * dx;
    rcost.lu = dt * R * du;
    rcost.ly = dt * S * dy;
    rcost.lxx = dt * Q;
    rcost.luu = dt * R;
    rcost.lux.setZero();
    rcost.lyy = dt * S;    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::terminal_cost(TCost&tcost, const State& x, int kend)
{
    const auto& dx = x - Xr->at(kend);
    tcost.Phi = dx.transpose() * Qf * dx;
    tcost.Phi *= 0.5;
}                                                   

template <typename T, size_t xs_, size_t us_, size_t ys_>
void QuadraticCost<T, xs_, us_, ys_>::terminal_cost_par(TCost&tcost, const State& x, int kend)
{  
    const auto& dx = x - Xr->at(kend);
    tcost.Phix = Qf * dx;
    tcost.Phixx = Qf;
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::running_cost(RCost& rcostdata, const State& x, const Contrl& u, 
                                                   const Output& y, T dt, int k)
{
    rcostdata.Zeros();
    for (auto cost_ptr:cost_ptrs)
    {
        rcostdata_temp.Zeros();
        cost_ptr->running_cost(rcostdata_temp, x, u, y, dt, k);     
        rcostdata.add(rcostdata_temp);
    }    
}                                                

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::running_cost_par(RCost& rcostdata, const State& x, const Contrl& u, 
                                                       const Output& y, T dt, int k)
{
    for (auto cost_ptr:cost_ptrs)
    {
        rcostdata_temp.Zeros();
        cost_ptr->running_cost_par(rcostdata_temp, x, u, y, dt, k);     
        rcostdata.add(rcostdata_temp);
    }    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::terminal_cost(TCost&tcostdata, const State& x, int kend)
{
    tcostdata.Zeros();
    for (auto cost_ptr:cost_ptrs)
    {
        tcostdata_temp.Zeros();
        cost_ptr->terminal_cost(tcostdata_temp, x, kend);
        tcostdata.add(tcostdata_temp);
    }    
}

template <typename T, size_t xs_, size_t us_, size_t ys_>
void CostContainer<T, xs_, us_, ys_>::terminal_cost_par(TCost&tcostdata, const State& x, int kend)
{
    for (auto cost_ptr:cost_ptrs)
    {
        tcostdata_temp.Zeros();
        cost_ptr->terminal_cost_par(tcostdata_temp, x, kend);
        tcostdata.add(tcostdata_temp);
    }    
}

// Explicit specilization for QuadraticCost when implemend separately in cpp
template class QuadraticCost<double,24,24,0>;
// template class QuadraticCost<float,24,24,0>;
template class CostContainer<double,24,24,0>;
