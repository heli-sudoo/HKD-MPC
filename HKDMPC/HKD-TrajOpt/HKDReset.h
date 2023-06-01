#ifndef HKDRESET_H
#define HKDRESET_H

#include "HSDDP_CPPTypes.h"
#include "casadiGen_HKD.h"
#include <cassert>
#include <vector>

using std::vector;
template <typename T>
class HKDReset
{
public:
    static const size_t xs = 24;
    static const size_t us = 24;
    static const size_t ys = 0;

private:
    VecM<T, 3> pos, eul, qleg, qleg_default, pf;
    VecM<T, 12> qdummy;
    MatMN<T, 3, 18> Jf; // foot Jacobian

public:
    HKDReset(/* args */);
    void resetmap(DVec<T> &xnext, DVec<T> &x, VecM<int, 4> &c, VecM<int, 4> &cn);
    void resetmap_partial(DMat<T> &Px, DVec<T> &x, VecM<int, 4> &c, VecM<int, 4> &cn);
};

template <typename T>
HKDReset<T>::HKDReset(/* args */)
{
    pos.setZero();
    eul.setZero();
    qleg.setZero();
    pf.setZero();
    qdummy.setZero();
    qleg_default << 0.0, -0.8, 1.7;
}

template <typename T>
void HKDReset<T>::resetmap(DVec<T> &xnext, DVec<T> &x, VecM<int, 4> &c, VecM<int, 4> &cn)
{
    assert((x.size() == 24));

    xnext.setZero(24);
    eul = x.head(3);
    pos = x.segment(3, 3);
    qdummy = x.tail(12);

    for (int l = 0; l < 4; l++)
    {
        // set to default joint angle if from stance to swing
        if (c[l] && (!cn[l]))
        {
            qdummy.segment(3 * l, 3) = qleg_default;
        }
        // compute foot placement if from swing to stance
        if ((!c[l]) && cn[l])
        {            
            qleg = qdummy.segment(3 * l, 3);
            pf.setZero();
            // compute foot position for foot foot_id
            T leg_id_T = static_cast<T>(l) + 1;
            vector<T *> arg_pos = {pos.data(), eul.data(), qleg.data(), &leg_id_T};
            vector<T *> res_pos = {pf.data()};
            casadi_interface(arg_pos, res_pos, pf.size(), compute_foot_position,
                             compute_foot_position_sparsity_out,
                             compute_foot_position_work);
            VecM<T,3> cmap;
            cmap << 1,1,0;
            qdummy.segment(3 * l, 3) = cmap.asDiagonal() * pf;
        }
    }
    xnext << x.head(12), qdummy;
}

template <typename T>
void HKDReset<T>::resetmap_partial(DMat<T> &Px, DVec<T> &x, VecM<int, 4> &c, VecM<int, 4> &cn)
{
    assert((x.size() == 24));
    Px.setIdentity(24, 24);

    eul = x.head(3);
    pos = x.segment(3, 3);
    qdummy = x.tail(12);

    for (int l = 0; l < 4; l++)
    {
        // set to zero mapping if from stance to swing
        if (c[l] && (!cn[l]))
        {
            Px.middleRows(12 + 3*l, 3).setZero();
        }       
        // use foot jacobian if from swing to stance         
        if ((!c[l]) && cn[l])
        {            
            qleg = qdummy.segment(3 * l, 3);
            Jf.setZero();            
            // compute foot Jacobian for foot foot_id
            vector<T *> arg_J = {pos.data(), eul.data(), qleg.data()};
            vector<T *> res_J = {Jf.data()};
            switch (l)
            {
            case 0: // FR
                casadi_interface(arg_J, res_J, Jf.size(), comp_foot_jacob_1,
                                 comp_foot_jacob_1_sparsity_out,
                                 comp_foot_jacob_1_work);
                break;
            case 1: // FL
                casadi_interface(arg_J, res_J, Jf.size(), comp_foot_jacob_2,
                                 comp_foot_jacob_2_sparsity_out,
                                 comp_foot_jacob_2_work);
                break;
            case 2: // HR
                casadi_interface(arg_J, res_J, Jf.size(), comp_foot_jacob_3,
                                 comp_foot_jacob_3_sparsity_out,
                                 comp_foot_jacob_3_work);
                break;
            case 3: // HL
                casadi_interface(arg_J, res_J, Jf.size(), comp_foot_jacob_4,
                                 comp_foot_jacob_4_sparsity_out,
                                 comp_foot_jacob_4_work);
                break;
            default:
                break;
            }            
            VecM<T,3> cmap;
            cmap << 1,1,0;
            Jf = cmap.asDiagonal()*Jf;
            // fill in Px using Jf appropriately
            Px.block(12+3*l, 0, 3, 3) = Jf.middleCols(3,3);
            Px.block(12+3*l, 3, 3, 3) = Jf.leftCols(3);            
            Px.block(12+3*l, 12,3, 12) = Jf.rightCols(12);
        }
    }
}
#endif
