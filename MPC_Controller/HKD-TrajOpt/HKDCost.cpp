#include "HKDCost.h"
#include "HSDDP_Utils.h"

template <typename T>
void HKDFootPlaceReg<T>::running_cost(RCost& rcost, const State& x, const Contrl& u, const Output& y, T dt, int k)
{
    pCoM_aug << x.segment(3,3), x.segment(3,3), x.segment(3,3), x.segment(3,3);
    prel = x.tail(12) - pCoM_aug;
    prel_r = Xr->at(k).tail(12);
    d_prel = prel - prel_r;
    rcost.l = .5 * d_prel.transpose() * Qfoot * d_prel;
    rcost.l *= dt;
}

template <typename T>
void HKDFootPlaceReg<T>::running_cost_par(RCost& rcost, const State& x, const Contrl& u, const Output& y, T dt, int k)
{
    pCoM_aug << x.segment(3,3), x.segment(3,3), x.segment(3,3), x.segment(3,3);
    prel = x.tail(12) - pCoM_aug;
    prel_r = Xr->at(k).tail(12);
    d_prel = prel - prel_r;
    rcost.lx = dt * dprel_dx.transpose() * Qfoot * d_prel;
    rcost.lxx = dt * dprel_dx.transpose() * Qfoot * dprel_dx;
}

template <typename T>
void HKDFootPlaceReg<T>::terminal_cost(TCost& tcost, const State& x, int kend)
{
    pCoM_aug << x.segment(3,3), x.segment(3,3), x.segment(3,3), x.segment(3,3);
    prel = x.tail(12) - pCoM_aug;
    prel_r = Xr->at(kend).tail(12);
    d_prel = prel - prel_r;
    tcost.Phi = 10 * d_prel.transpose() * Qfoot * d_prel;
    
    // turning
    // tcost.Phi = 5 * d_prel.transpose() * Qfoot * d_prel;
}

template <typename T>
void HKDFootPlaceReg<T>::terminal_cost_par(TCost& tcost, const State& x, int kend)
{
    pCoM_aug << x.segment(3,3), x.segment(3,3), x.segment(3,3), x.segment(3,3);
    prel = x.tail(12) - pCoM_aug;
    prel_r = Xr->at(kend).tail(12);
    d_prel = prel - prel_r;
    tcost.Phix = 20 * dprel_dx.transpose() * Qfoot * d_prel;
    tcost.Phixx = 20 * dprel_dx.transpose() * Qfoot * dprel_dx;

    // turning
    // tcost.Phix = 10 * dprel_dx.transpose() * Qfoot * d_prel;
    // tcost.Phixx = 10 * dprel_dx.transpose() * Qfoot * dprel_dx;
}

template <typename T>
void PrevSolution_Reg<T>::running_cost(RCost&rcost, const State& x, const Contrl& u, const Output& y, T dt, int k)
{
    unused_ignore(y);
    unused_ignore(dt);
    const auto& dx = x - Xr->at(k);
    const auto& du = u - Ur->at(k);
    rcost.l = 0.5*dx.transpose() * Q * dx;
    rcost.l += 0.5*du.transpose() * R * du;      
}

template <typename T>
void PrevSolution_Reg<T>::running_cost_par(RCost&rcost, const State& x, const Contrl& u, const Output& y, T dt, int k)
{
    unused_ignore(y);
    unused_ignore(dt);
    const auto& dx = x - Xr->at(k);
    const auto& du = u - Ur->at(k);
    rcost.lx = Q * dx;
    rcost.lu = R * du;
    rcost.lxx = Q;
    rcost.luu = R;     
}

template <typename T>
void PrevSolution_Reg<T>::terminal_cost(TCost&tcost, const State& x, int kend)
{
    const auto& dx = x - Xr->at(kend);
    tcost.Phi = dx.transpose() * Qf * dx;
    tcost.Phi *= 0.5; 
}

template <typename T>
void PrevSolution_Reg<T>::terminal_cost_par(TCost&tcost, const State& x, int kend)
{
    const auto& dx = x - Xr->at(kend);
    tcost.Phix = Qf * dx;
    tcost.Phixx = Qf;
   
}

template class HKDFootPlaceReg<double>;
template class PrevSolution_Reg<double>;

