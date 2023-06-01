#ifndef HKDCOST_H
#define HKDCOST_H
#include "SinglePhaseInterface.h"
#include "QuadReference.h"
#include "HKDModel.h"

template <typename T>
class HKDTrackingCost : public QuadraticTrackingCost<T, HKD::xs, HKD::us, HKD::ys>
{
public:
    HKDTrackingCost(VecM<int, 4> contact) : QuadraticTrackingCost<T, HKD::xs, HKD::us, HKD::ys>()
    {
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_eul(1, 4, 5);
        VecM<T, 3> q_pos(1, 1, 30);
        VecM<T, 3> q_omega(.2, .2, .2);
        VecM<T, 3> q_v(4, 1, .5);
        VecM<T, 12> q_qJ;                        
        q_qJ << VecM<T, 3>::Constant(.2 * (1 - contact[0])),
                VecM<T, 3>::Constant(.2 * (1 - contact[1])),
                VecM<T, 3>::Constant(.2 * (1 - contact[2])),
                VecM<T, 3>::Constant(.2 * (1 - contact[3]));
       this->Q.diagonal() << q_eul, q_pos, q_omega, q_v, q_qJ;
        
        /* Terminal state weighting matrix */
        VecM<T, 24> scale;
        scale << 1, 1, 2, 1, 1, 20, 
                .3, .3, .3,  1, 3, 1, 
                .01 * VecM<T, 12>::Ones();
        this->Qf = 20 * scale.asDiagonal() * this->Q;

        /* Control weighting matrices */
        VecM<T, 24> r;        
        r.template head<12>() = .2 * VecM<T,12>::Ones();      // GRF
        r.template tail<12>() = .1 * VecM<T,12>::Ones();  // Commanded joint vel       
        this->R = r.asDiagonal();      
    }       
};

template <typename T>
class HKDFootPlaceReg : public CostBase<T, HKD::xs, HKD::us, HKD::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,HKD::xs, HKD::us, HKD::ys>::State;
    using typename CostBase<T,HKD::xs, HKD::us, HKD::ys>::Contrl;
    using typename CostBase<T,HKD::xs, HKD::us, HKD::ys>::Output;
    using typename CostBase<T,HKD::xs, HKD::us, HKD::ys>::RCost;
    using typename CostBase<T,HKD::xs, HKD::us, HKD::ys>::TCost;    

    HKDFootPlaceReg(const VecM<int, 4>& contact) : 
    CostBase<T, HKD::xs, HKD::us, HKD::ys>("Foot regularization"){
        Qfoot.setZero();
        dprel_dx.setZero();
        Qfoot.diagonal() << 3*contact[0],contact[0],0,
                            3*contact[1],contact[1],0,
                            3*contact[2],contact[2],0,
                            3*contact[3],contact[3],0;              
        dprel_dx.middleCols(3,3) << -contact[0]*Mat3<T>::Identity(),
                                    -contact[1]*Mat3<T>::Identity(),
                                    -contact[2]*Mat3<T>::Identity(),
                                    -contact[3]*Mat3<T>::Identity();
        dprel_dx.block(0, 12, 3, 3) << contact[0]*Mat3<T>::Identity();
        dprel_dx.block(3, 15, 3, 3) << contact[1]*Mat3<T>::Identity();
        dprel_dx.block(6, 18, 3, 3) << contact[2]*Mat3<T>::Identity();
        dprel_dx.block(9, 21, 3, 3) << contact[3]*Mat3<T>::Identity(); 
        
        Qfoot *= 20;

        quad_reference = nullptr;
        quad_astate = nullptr;
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void terminal_cost(TCost&, const State& x, float tend=0) override;

    void terminal_cost_par(TCost&, const State&x, float tend=0) override;

    void set_quad_reference(QuadReference* quad_reference_in){
        quad_reference = quad_reference_in;                        
    }

private:
    MatMN<T, 12, 12> Qfoot;
    MatMN<T, 12, 24> dprel_dx;

    VecM<double, 12> prel_r;
    VecM<double, 3> pCoM_r;
    VecM<T, 12> prel;
    VecM<T, 3> pCoM;    
    VecM<T, 12> d_prel;

    QuadReference* quad_reference;
    QuadAugmentedState* quad_astate;
};

#endif