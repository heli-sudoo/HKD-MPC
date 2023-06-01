#ifndef HKDMODEL_H
#define HKDMODEL_H

#include "HSDDP_CPPTypes.h"
#include <vector>
#include <iostream>
#include "casadiGen_HKD.h"
#include "HSDDP_Utils.h"

namespace HKD
{
    const size_t xs = 24;
    const size_t us = 24;
    const size_t ys = 0;

    template <typename T>
    class Model
    {
    public:
        typedef VecM<T, xs> StateType;
        typedef VecM<T, us> ContrlType;
        typedef VecM<T, ys> OutputType;
        typedef VecM<T, 12> JointType;

        typedef MatMN<T, xs, xs> StateMap;
        typedef MatMN<T, xs, us> ContrlMap;
        typedef MatMN<T, ys, xs> OutputMap;
        typedef MatMN<T, ys, us> DirectMap;

        typedef VecM<int, 4> CtactStatusType;

    public:
        void dynamics(StateType &xnext, OutputType &y,
                      StateType &x, ContrlType &u, T t,
                      CtactStatusType &ctact_status, T &dt)
        {
            (void)(y);
            (void)(t);
            VecM<T, 4> ctact_status_T = ctact_status.cast<T>();
            std::vector<T *> arg = {x.data(), u.data(), &dt, ctact_status_T.data()};
            std::vector<T *> res = {xnext.data()};
            casadi_interface(arg, res, xnext.size(), hkinodyn,
                             hkinodyn_sparsity_out,
                             hkinodyn_work);
        }
        void dynamics_partial(StateMap &A, ContrlMap &B, OutputMap &C, DirectMap &D,
                              StateType &x, ContrlType &u, T t,
                              CtactStatusType &ctact_status, T &dt)
        {
            (void)(C);
            (void)(D);
            (void)(t);
            B.setZero();
            A.setZero();
            VecM<T, 4> ctact_status_T = ctact_status.cast<T>();
            std::vector<T *> arg = {x.data(), u.data(), &dt, ctact_status_T.data()};
            std::vector<T *> res = {A.data(), B.data()};
            casadi_interface(arg, res, A.size(), hkinodyn_par,
                             hkinodyn_par_sparsity_out,
                             hkinodyn_par_work);
        }
    };
}

template <typename D1, typename D2>
void compute_hkd_state(D1 &eul, D1 &pos,
                       D2 &qJ,
                       D2 &qdummy, const VecM<int, 4> &c)
{
    typedef typename D1::Scalar T;
    // initialize qdummy based on whether it is joint angle or foot position
    for (int l = 0; l < 4; l++)
    {
        // get the joint angle for leg l
        D1 qleg = qJ.template segment(3 * l, 3);

        // if swing, set qdummy to joint angle
        if (c[l] == 0)
        {
            qdummy.segment(3 * l, 3) = qleg;
        }
        // if stance, set qdummy to foot location
        else
        {
            T foot_id_T = static_cast<T>(l) + 1;
            // compute foot position for foot foot_id
            VecM<T, 3> pf_T;
            std::vector<T *> arg_pos = {pos.data(), eul.data(), qleg.data(), &foot_id_T};
            std::vector<T *> res_pos = {pf_T.data()};
            casadi_interface(arg_pos, res_pos, pf_T.size(), compute_foot_position,
                             compute_foot_position_sparsity_out,
                             compute_foot_position_work);
            qdummy.segment(3 * l, 3) = pf_T;
        }
    }
}

#endif