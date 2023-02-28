#ifndef HKDCOST_H
#define HKDCOST_H
#include "CostBase.h"
#include <deque>

using std::deque;

template <typename T>
class HKDTrackingCost : public QuadraticCost<T, 24, 24, 0>
{
public:
    HKDTrackingCost(VecM<int, 4> contact) : QuadraticCost<T, 24, 24, 0>()
    {
        Qeul.setZero();
        Qpos.setZero();
        Qw.setZero();
        Qv.setZero();
        QqJ.setZero();
                
        Qeul.diagonal() << 1, 4, 5;
        Qpos.diagonal() << .2, .2, 40;
        Qw.diagonal() << .2, .2, .2;
        Qv.diagonal() << 4, 1, .1;
        QqJ.diagonal() << VecM<T, 3>::Constant(.2 * (1 - contact[0])),
                        VecM<T, 3>::Constant(.2 * (1 - contact[1])),
                        VecM<T, 3>::Constant(.2 * (1 - contact[2])),
                        VecM<T, 3>::Constant(.2 * (1 - contact[3]));
        this->Q.topLeftCorner(3, 3) << Qeul;
        this->Q.block(3, 3, 3, 3) << Qpos;
        this->Q.block(6, 6, 3, 3) << Qw;
        this->Q.block(9, 9, 3, 3) << Qv;
        this->Q.bottomRightCorner(12, 12) << QqJ;

        /* Terminal state weighting matrix */
        VecM<T, 24> scale;
        scale << 1, 1, 2, 
                  .2, .2, 20, 
                  .3, .3, .3, 
                  1, 3, 1, 
                  .01 * VecM<T, 12>::Ones();
        this->Qf = 20 * scale.asDiagonal() * this->Q;

        /* Control weighting matrices */
        this->R.setZero();
        this->R.topLeftCorner(12, 12) = .2 * MatMN<T, 12, 12>::Identity();      // GRF
        this->R.bottomRightCorner(12, 12) = .01 * MatMN<T, 12, 12>::Identity();  // Commanded joint vel
    }

private:
    MatMN<T, 3, 3> Qeul;
    MatMN<T, 3, 3> Qpos;
    MatMN<T, 3, 3> Qw;
    MatMN<T, 3, 3> Qv;
    MatMN<T, 12, 12> QqJ;
    
};

template <typename T>
class HKDFootPlaceReg : public CostBase<T, 24, 24, 0>
{
public:
    using typename CostBase<T,24,24,0>::State;
    using typename CostBase<T,24,24,0>::Contrl;
    using typename CostBase<T,24,24,0>::Output;
    using typename CostBase<T,24,24,0>::RCost;
    using typename CostBase<T,24,24,0>::TCost;

    HKDFootPlaceReg(VecM<int, 4> contact) : CostBase<T, 24, 24, 0>(){
        Qfoot.setZero();
        dprel_dx.setZero();
        Qfoot.diagonal() << 3*contact[0],contact[0],0,
                            3*contact[1],contact[1],0,
                            3*contact[2],contact[2],0,
                            3*contact[3],contact[3],0;    
        // Qfoot.diagonal() << 0,contact[0],0,
        //                     0,contact[1],0,
        //                     0,contact[2],0,
        //                     0,contact[3],0;       
        dprel_dx.middleCols(3,3) << -contact[0]*Mat3<T>::Identity(),
                                    -contact[1]*Mat3<T>::Identity(),
                                    -contact[2]*Mat3<T>::Identity(),
                                    -contact[3]*Mat3<T>::Identity();
        dprel_dx.block(0, 12, 3, 3) << contact[0]*Mat3<T>::Identity();
        dprel_dx.block(3, 15, 3, 3) << contact[1]*Mat3<T>::Identity();
        dprel_dx.block(6, 18, 3, 3) << contact[2]*Mat3<T>::Identity();
        dprel_dx.block(9, 21, 3, 3) << contact[3]*Mat3<T>::Identity(); 
        
        Qfoot *= 20;
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void terminal_cost(TCost&, const State& x, int kend=0) override;

    void terminal_cost_par(TCost&, const State&x, int kend=0) override;

    void set_reference(deque<VecM<T, 24>>* Xr_in){
        Xr = Xr_in;
    }

private:
    MatMN<T, 12, 12> Qfoot;
    MatMN<T, 12, 24> dprel_dx;

    VecM<T, 12> prel_r;
    VecM<T, 12> prel;
    VecM<T, 12> pCoM_aug;
    VecM<T, 12> d_prel;
    deque<VecM<T, 24>>* Xr;
};

template<typename T>
class PrevSolution_Reg : public CostBase<T,24,24,0>
{
public:
    using typename CostBase<T,24,24,0>::State;
    using typename CostBase<T,24,24,0>::Contrl;
    using typename CostBase<T,24,24,0>::Output;
    using typename CostBase<T,24,24,0>::RCost;
    using typename CostBase<T,24,24,0>::TCost;

    PrevSolution_Reg(VecM<int, 4> contact) : CostBase<T, 24, 24, 0>()
    {
        Qfoot.setZero();
        Qfoot.diagonal() << contact[0],contact[0],0,
                            contact[1],contact[1],0,
                            contact[2],contact[2],0,
                            contact[3],contact[3],0;
        Q.setIdentity();
        Q.bottomRightCorner(12,12) = 10*Qfoot;
        Q *= 0.1;
        Qf = 5 * Q;
        R.setZero();
    }

public:
    
    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, int k=0) override;

    void terminal_cost(TCost&, const State& x, int kend=0) override;

    void terminal_cost_par(TCost&, const State&x, int kend=0) override;

    void set_state_reference(deque<VecM<T, 24>>* Xr_in){
        Xr = Xr_in;
    }

    void set_control_reference(deque<VecM<T, 24>>* Ur_in)
    {
        Ur = Ur_in;
    }

    void set_maximum_step(int max_step_in){max_step = max_step_in;}

private:
    deque<VecM<T, 24>>* Xr = nullptr;
    deque<VecM<T, 24>>* Ur = nullptr;

    MatMN<T,24,24> Q;
    MatMN<T,24,24> R;
    MatMN<T,24,24> Qf;
    MatMN<T,12,12> Qfoot;
    int max_step;
};
#endif