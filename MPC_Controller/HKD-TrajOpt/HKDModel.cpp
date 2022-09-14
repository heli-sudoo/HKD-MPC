// #include "HKDModel.h"
// #include "CasadiGen.h"
// #include "HSDDP_Utils.h"
// template<typename T>
// void HKDModel<T>::dynamics(StateType& xnext, OutputType& y,
//                            StateType& x, ContrlType& u, 
//                            CtactStatusType& ctact_status, T& dt)
// {
//     unused_ignore(y);
//     update_foot_loc(x, ctact_status);
//     VecM<T, 4> ctact_status_T = ctact_status.cast<T>();
//     vector<T *> arg = {x.data(), u.data(), &dt, foot_loc_world.data(), ctact_status_T.data()};
//     vector<T *> res = {xnext.data()};
//     casadi_interface(arg, res, xnext.size(), hkinodyn, 
//                      hkinodyn_sparsity_out,
//                      hkinodyn_work);
// }

// template<typename T>
// void HKDModel<T>::dynamics_partial(StateMap& A, ContrlMap& B, OutputMap& C, DirectMap& D,
//                                    StateType& x, ContrlType& u, 
//                                    CtactStatusType& ctact_status, T& dt)
// {
//     unused_ignore(C);
//     unused_ignore(D);
//     B.setZero();
//     A.setZero();
//     update_foot_loc(x, ctact_status);
//     VecM<T, 4> ctact_status_T = ctact_status.cast<T>();
//     vector<T *> arg = {x.data(), u.data(), &dt, foot_loc_world.data(), ctact_status_T.data()};
//     vector<T *> res = {A.data(), B.data()};
//     casadi_interface(arg, res, A.size(), hkinodyn_par,
//                      hkinodyn_par_sparsity_out,
//                      hkinodyn_par_work);
// }   

// template<typename T>
// void HKDModel<T>::update_foot_loc(StateType&x, CtactStatusType& ctact_status)
// {
//     /* eul and pos are flipped in the computation of foot position and Jacobian */
//     VecM<T, 18> q;
//     q.head(3) = x.segment(3,3);
//     q.segment(3,3) = x.head(3);
//     q.tail(12) = x.tail(12);
//     /* compute foot positions */
//     vector<T *> arg_pos = {q.data()};
//     vector<T *> res_pos = {foot_loc_world_temp.data()};
//     casadi_interface(arg_pos, res_pos, foot_loc_world_temp.size(), compute_foot_positions, 
//                 compute_foot_positions_sparsity_out, 
//                 compute_foot_positions_work);
//     /* set foot height to be always zero */
//     foot_loc_world_temp(2) = 0;
//     foot_loc_world_temp(5) = 0;
//     foot_loc_world_temp(8) = 0;
//     foot_loc_world_temp(11) = 0;

//     /* update foot position for swing foot */
//     for (size_t foot = 0; foot < 4; foot++)
//     {
//         if (!ctact_status[foot])
//         {
//             foot_loc_world.segment(3*foot, 3) = foot_loc_world_temp.segment(3*foot, 3);
//         }
        
//     }        
// }

// template class HKDModel<double>;