#include "HKDConstraints.h"
#include "HSDDP_Utils.h"
#include "casadiGen_HKD.h"
#include <cassert>

template <typename T>
GRFConstraint<T>::GRFConstraint(const VecM<int, 4> &ctact)
    : PathConstraintBase<T, 24, 24, 0>("GRF")
{
    ctact_status = ctact;
    int n_constrained_foot = ctact.cwiseEqual(1).cast<int>().sum();
    this->update_constraint_size(5 * n_constrained_foot);
   
    auto &mu = mu_fric;
    A.setZero(this->size, 24);
    b.setZero(this->size);
    MatMN<T, 5, 3> A_leg;
    A_leg << 0, 0, 1,
        -1, 0, sqrt(mu),
        1, 0, sqrt(mu),
        0, -1, sqrt(mu),
        0, 1, sqrt(mu);

    int i(0);
    for (size_t leg = 0; leg < 4; ++leg)
    {
        if (ctact_status[leg] > 0)
        {
            A.block(5 * i, 3 * leg, 5, 3) = A_leg;
            ++i;
        }               
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

    (void)(x);
    (void)(y);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {
        this->data[k][i].g = A.row(i) * u;        
    }
    this->update_max_violation(k);
}

template <typename T>
void GRFConstraint<T>::compute_partial(const State &x, const Contrl &u, const Output &y, int k)
{
    (void)(x);
    (void)(y);
    (void)(u);

    for (size_t i = 0; i < this->data[k].size(); i++)
    {        
        this->data[k][i].gu = A.row(i).transpose();
    }
}

template <typename T>
TouchDownConstraint<T>::TouchDownConstraint(const VecM<int, 4> &impact_status_in)
    : TerminalConstraintBase<T, 24>("TouchDwon")
{
    ground_height = 0;
    impact_status = impact_status_in;
    impact_foot_ids = find_eigen(impact_status, 1);
    this->update_constraint_size(impact_foot_ids.size());
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
    eul = x.template head<3>();
    pos = x.template segment<3>(3);
    qdummy = x.template tail<12>();

    /* iterate through all constraints */
    for (size_t i = 0; i < this->size; i++)
    {
        int foot_id = impact_foot_ids[i];
        T foot_id_T = static_cast<T>(foot_id)+1;

        qleg = qdummy.template segment<3>(3 * foot_id);
        std::vector<T *> arg_pos = {pos.data(), eul.data(), qleg.data(), &foot_id_T};
        std::vector<T *> res_pos = {pFoot.data()};

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
    eul = x.template head<3>();
    pos = x.template segment<3>(3);
    qdummy = x.template tail<12>();

    // iterate through all constraints
    for (size_t i = 0; i < this->size; i++)
    {
        int foot_id = impact_foot_ids[i];

        qleg = qdummy.template segment<3>(3 * foot_id);
    
        // compute foot Jacobian
        Jf.setZero();
        std::vector<T *> arg_J = {pos.data(), eul.data(), qleg.data()};
        std::vector<T *> res_J = {Jf.data()};
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
        this->data[i].hx.template head<3>() = Jz.template segment<3>(3);
        this->data[i].hx.template segment<3>(3) = Jz.template head<3>();
        this->data[i].hx.template tail<12>() = Jz.template tail<12>();
    }
}

template class GRFConstraint<double>;
template class TouchDownConstraint<double>;