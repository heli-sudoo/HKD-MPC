#include "HKDReference.h"

/* 
    @brief: Get a quadruped reference state and control at time t, and converts to the hkd state refrence
    @params: t: relative time w.r.t. the beginning of the current reference (starting time is always zero)
*/

void HKDSinglePhaseReference::get_reference_at_t(VecM<double, 24>& xt, VecM<double, 24>& ut, VecM<double, 0>& yt, float t)
{
    
    (void) (yt);
    
    get_reference_at_t(xt, t);

    ut.head<12>() = quad_state_t_ptr->grf;       // grf reference
    ut.tail<12>() = quad_state_t_ptr->qJd;       // commanded joint vel reference    
}

/* 
    @brief: Get a quadruped reference state at time t, and converts to the hkd state refrence
    @params: t: relative time w.r.t. the beginning of the current reference (starting time is always zero)
*/

void HKDSinglePhaseReference::get_reference_at_t(VecM<double, 24>& xt, float t)
{
    // Check nullptr
    if (nullptr == quad_ref_)
    {
        printf("error: quad_ref is nullptr \n");
        printf("return \n");
        return;
    }

    quad_state_t_ptr = quad_ref_->get_a_reference_ptr_at_t(t);

    if (nullptr == quad_state_t_ptr)
    {
        printf("error: quad_state_ref_ptr is nullptr \n");
        return;
    }

    xt.head<12>() = quad_state_t_ptr->body_state;   //body state: euler, pos, ang, vel

    for (int leg = 0; leg < 4; leg++)
    {
        // If in stance
        if (quad_state_t_ptr->contact[leg] > 0)
        {
            xt.segment<3>(12+3*leg) = quad_state_t_ptr->foot_placements.segment<3>(3*leg);
        }
        // If in swing
        else
        {
            xt.segment<3>(12+3*leg) = quad_state_t_ptr->qJ.segment<3>(3*leg);
        }                
    }
}
