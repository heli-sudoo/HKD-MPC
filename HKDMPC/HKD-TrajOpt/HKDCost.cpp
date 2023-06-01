#include "HKDCost.h"
#include "HSDDP_Utils.h"

template <typename T>
void HKDFootPlaceReg<T>::running_cost(RCost& rcost, const State& x, const Contrl& u, const Output& y, T dt, float t)
{
    (void) (y);
    (void) (u);  
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    pCoM << x.template segment<3>(3);
    pCoM_r << quad_astate->body_state.segment<3>(3);
                 
    prel = x.tail(12) - pCoM.template replicate<4,1>();
    prel_r = quad_astate->foot_placements - pCoM_r.replicate<4,1>();
    d_prel = prel - prel_r.cast<T>();
    rcost.l = .5 * d_prel.transpose() * Qfoot * d_prel;
    rcost.l *= dt;    
}

template <typename T>
void HKDFootPlaceReg<T>::running_cost_par(RCost& rcost, const State& x, const Contrl& u, const Output& y, T dt, float t)
{
    (void) (y);
    (void) (u);
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    pCoM << x.template segment<3>(3);
    pCoM_r << quad_astate->body_state.segment<3>(3);

    prel = x.tail(12) - pCoM.template replicate<4,1>();
    prel_r = quad_astate->foot_placements - pCoM_r.replicate<4,1>();
    d_prel = prel - prel_r.cast<T>();
    rcost.lx = dt * dprel_dx.transpose() * Qfoot * d_prel;
    rcost.lxx = dt * dprel_dx.transpose() * Qfoot * dprel_dx;
}

template <typename T>
void HKDFootPlaceReg<T>::terminal_cost(TCost& tcost, const State& x, float tend)
{
    quad_astate = quad_reference->get_a_reference_ptr_at_t(tend);

    pCoM << x.template segment<3>(3);
    pCoM_r << quad_astate->body_state.segment<3>(3);

    prel = x.tail(12) - pCoM.template replicate<4,1>();
    prel_r = quad_astate->foot_placements - pCoM_r.replicate<4,1>();
    d_prel = prel - prel_r.cast<T>();
    tcost.Phi = 10 * d_prel.transpose() * Qfoot * d_prel;  
}

template <typename T>
void HKDFootPlaceReg<T>::terminal_cost_par(TCost& tcost, const State& x, float tend)
{
    quad_astate = quad_reference->get_a_reference_ptr_at_t(tend);
    
    pCoM << x.template segment<3>(3);
    pCoM_r << quad_astate->body_state.segment<3>(3);

    prel = x.tail(12) - pCoM.template replicate<4,1>();
    prel_r = quad_astate->foot_placements - pCoM_r.replicate<4,1>();
    d_prel = prel - prel_r.cast<T>();
    tcost.Phix = 20 * dprel_dx.transpose() * Qfoot * d_prel;
    tcost.Phixx = 20 * dprel_dx.transpose() * Qfoot * dprel_dx;
  
}


template class HKDFootPlaceReg<double>;

