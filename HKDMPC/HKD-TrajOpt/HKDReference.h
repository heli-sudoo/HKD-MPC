#ifndef HKDREFERENCE_H
#define HKDREFERENCE_H

#include <deque>
#include <memory>
#include <algorithm>
#include <cassert>
#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"
#include "SinglePhaseInterface.h"
#include "QuadReference.h"


        
class HKDSinglePhaseReference : public SinglePhaseReferenceAbstract<24,24,0>
{
public:
    HKDSinglePhaseReference() {}

    void get_reference_at_t(VecM<double, 24>& xt,  float t) override;

    void get_reference_at_t(VecM<double, 24>& xt, VecM<double, 24>& ut, VecM<double, 0>& yt, float t) override;

    void set_quadruped_reference(QuadReference* quad_ref) {quad_ref_ = quad_ref;}


private:
    QuadReference* quad_ref_ = nullptr;  

    QuadAugmentedState* quad_state_t_ptr = nullptr;  

};



#endif //HKDREFERENCE_H