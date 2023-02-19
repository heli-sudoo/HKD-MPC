#include "HKDConstraints.h"
#include "HSDDP_Utils.h"
#include "CasadiGen.h"
#include <cassert>

template <typename T>
GRFConstraint<T>::GRFConstraint(const VecM<int, 4> &ctact)
    : PathConstraintBase<T, 24, 24, 0>("GRF")
{
    update_contact_status(ctact);
    auto &mu = mu_fric;
    A.setZero(this->size, 24);
    b.setZero(this->size);
    MatMN<T, 5, 3> A_leg;
    A_leg << 0, 0, 1,
        -1, 0, sqrt(mu),
        1, 0, sqrt(mu),
        0, -1, sqrt(mu),
        0, 1, sqrt(mu);

    for (int i = 0; i < ctact_foot_ids.size(); i++)
    {
        int leg = ctact_foot_ids[i];
        A.block(5 * i, 3 * leg, 5, 3) = A_leg;
    }
}

template <typename T>
void GRFConstraint<T>::compute_violation(const State &x, const Contrl &u, const Output &y, int k)
{
    if (this->data.empty())
    {
        printf("Constraint data is empty \n");
        printf("Creating constaint data \n");
        this->create_data();
    }

    unused_ignore(x);
    unused_ignore(y);

    for (int i = 0; i < this->data[k].size(); i++)
    {
        this->data[k][i].g = A.row(i) * u;        
    }
    this->update_max_violation(k);
}

template <typename T>
void GRFConstraint<T>::compute_partial(const State &x, const Contrl &u, const Output &y, int k)
{
    unused_ignore(x);
    unused_ignore(y);

    for (int i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gu = A.row(i).transpose();
    }
}

template <typename T>
TouchDownConstraint<T>::TouchDownConstraint(VecM<int, 4> &impact_status_)
    : TerminalConstraintBase<T, 24>("TouchDwon")
{
    ground_height = 0;
    update_impact_status(impact_status_);
}

template <typename T>
void TouchDownConstraint<T>::compute_violation(const State &x)
{
    if (this->size != this->data.size())
    {
        printf("constraint size and constraint data size do not match \n");
        return;
    }

    // If constraint data is empty, create data
    if (this->size == 0)
    {
        printf("constraint size is zero \n");
        return;
    }

    /* compute foot positions using the casadi interface */
    VecM<T, 3> pFoot;
    VecM<T, 3> eul, pos, qleg;
    VecM<T, 12> qdummy;

    eul = x.head(3);
    pos = x.segment(3, 3);
    qdummy = x.tail(12);

    // iterate through all constraints
    for (int i = 0; i < this->size; i++)
    {
        int foot_id = impact_foot_ids[i];
        T foot_id_T = static_cast<T>(foot_id)+1;

        qleg = qdummy.segment(3 * foot_id, 3);
        vector<T *> arg_pos = {pos.data(), eul.data(), qleg.data(), &foot_id_T};
        vector<T *> res_pos = {pFoot.data()};

        // compute foot position
        casadi_interface(arg_pos, res_pos, pFoot.size(), compute_foot_position,
                         compute_foot_position_sparsity_out,
                         compute_foot_position_work);
        
        // compute constraint violation
        this->data[i].h = pFoot[2] - ground_height;
       
    }
    this->update_max_violation();
}

template <typename T>
void TouchDownConstraint<T>::compute_partial(const State &x)
{
    /* compute foot positions using the casadi interface */    
    VecM<T, 3> eul, pos, qleg;
    VecM<T, 12> qdummy;
    MatMN<T, 3, 18> Jf;
    VecM<T, 18> Jz;

    eul = x.head(3);
    pos = x.segment(3, 3);
    qdummy = x.tail(12);

    // iterate through all constraints
    for (int i = 0; i < this->size; i++)
    {
        int foot_id = impact_foot_ids[i];
        T foot_id_T = static_cast<T>(foot_id)+1;

        qleg = qdummy.segment(3 * foot_id, 3);
    
        // compute foot Jacobian
        Jf.setZero();
        vector<T *> arg_J = {pos.data(), eul.data(), qleg.data()};
        vector<T *> res_J = {Jf.data()};
        switch (foot_id)
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
        Jz = Jf.bottomRows(1).transpose();
        this->data[i].hx.head(3) = Jz.segment(3, 3);
        this->data[i].hx.segment(3, 3) = Jz.head(3);
        this->data[i].hx.tail(12) = Jz.tail(12);
    }
}

template class GRFConstraint<double>;
template class TouchDownConstraint<double>;