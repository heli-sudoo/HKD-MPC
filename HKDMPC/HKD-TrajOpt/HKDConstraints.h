#ifndef HKDCONSTRAINT_H
#define HKDCONSTRAINT_H

#include <vector>
#include "HSDDP_CPPTypes.h"
#include "ConstraintsBase.h"
#include "HSDDP_Utils.h"

#include "mip_hopping.h"

#ifdef MIP_HOPPING
#include "utilities.h"
#endif

template<typename T>
class GRFConstraint:public PathConstraintBase<T,24,24,0>
{
private:
    using typename PathConstraintBase<T,24,24,0>::State;
    using typename PathConstraintBase<T,24,24,0>::Contrl;
    using typename PathConstraintBase<T,24,24,0>::Output;

    T mu_fric = 0.7;
    DMat<T> A;
    DVec<T> b;

    #ifdef MIP_HOPPING
    VecM<T, 3> center_point;
    VecM<T, 3> eul_terrain;
    #endif

public:
    GRFConstraint(const VecM<int, 4>& ctact_);

    void set_friction_coefficient(T mu_fric_in){
        mu_fric = mu_fric_in;
    }
      
    void compute_violation(const State&, const Contrl&, const Output&, int k) override;

    void compute_partial(const State&, const Contrl&, const Output&, int k) override;

    void set_terrain_info(const Vec3<T> center_point_in, const Vec3<T> eul_terrain_in){
        center_point = center_point_in;
        eul_terrain = eul_terrain_in;
    }

public:
    VecM<int, 4> ctact_status;
    std::vector<int> ctact_foot_ids;
};

template<typename T>
class TouchDownConstraint:public TerminalConstraintBase<T,24>
{
private:
    using typename TerminalConstraintBase<T,24>::State;    

public:
    TouchDownConstraint(const VecM<int, 4> &);

    #ifdef MIP_HOPPING
    void update_ground_height();
    #else
    void update_ground_height(T gheight_in){
        ground_height = gheight_in;
    }
    #endif   

    void compute_violation(const State&);

    void compute_partial(const State&);

    void set_terrain_info(const Vec3<T> center_point_in, const Vec3<T> eul_terrain_in){
        center_point = center_point_in;
        eul_terrain = eul_terrain_in;
        Mat3<T> R;
        EulerZYX_2_SO3(eul_terrain, R);
        plane_coefficients = R * Vec3<T>(0, 0, 1);
    }

private:
    VecM<int, 4> impact_status;
    T ground_height;
    std::vector<int>impact_foot_ids;

    VecM<T, 3> pFoot;
    VecM<T, 3> eul, pos, qleg;
    VecM<T, 12> qdummy;
    MatMN<T, 3, 18> Jf;
    VecM<T, 18> Jz;

    #ifdef MIP_HOPPING
    Vec3<T> center_point;
    Vec3<T> eul_terrain;
    Vec3<T> plane_coefficients;
    #endif

};

template<typename T>
class SwingConstraint:public PathConstraintBase<T, 24, 24, 0>
{
private:
    using typename PathConstraintBase<T,24,24,0>::State;
    using typename PathConstraintBase<T,24,24,0>::Contrl;
    using typename PathConstraintBase<T,24,24,0>::Output;

public:
    SwingConstraint(VecM<int, 4>& swing_status_):
        PathConstraintBase<T, 24, 24, 0>("SwingConstraint"){
        ground_height = 0;
        update_swing_status(swing_status_);
        foot_positions.setZero();
        Js_in_z.setZero();
        J_trans.setZero();
    }

    void update_ground_height(T gheight_){
        ground_height = gheight_;
    }

    void update_swing_status(VecM<int, 4> &swing_status_){
        swing_foot_ids.clear();
        swing_status = swing_status_;
        swing_foot_ids = find_eigen(swing_status, 1);
        this->update_constraint_size(swing_foot_ids.size());
    }
    void compute_violation(const State&, const Contrl&, const Output&, int k) override;

public:
    VecM<int, 4> swing_status;
    std::vector<int> swing_foot_ids;
    T ground_height;
    VecM<T, 12> foot_positions;
    MatMN<T, 4, 18> Js_in_z;    
    VecM<T, 18> J_trans;
};


#endif // HKDCONSTRAINT_H