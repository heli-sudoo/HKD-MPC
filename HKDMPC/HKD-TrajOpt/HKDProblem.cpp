#include <memory>
#include "HKDProblem.h"
#include "HKDCost.h"
#include "HSDDP_Utils.h"
#include "HKDConstraints.h"
#include "TrajectoryManagement.h"
#include "SinglePhase.h"
#include <functional>            // std::bind, and std::placeholder
#include <tabulate/table.hpp>   // pretty print of tables
#include "utilities.h"

namespace pc = std::placeholders; // predefined position placeholder for bind

template <typename T>
void HKDProblem<T>::initialization()
{
    printf("Initializing HKDProblem ... \n\n");

    /* Initialize the reference */
    quad_ref_ptr->initialize(plan_duration);

    /* Initialize the hkd_reference interface */
    hkd_reference.set_quadruped_reference(quad_ref_ptr);

    /* Initialize the phase timing and contact information */    
    VecM<int, 4> contact_prev, contact_cur;
    VecM<double, 4> contact_duration;    
    contact_prev.setZero();
    contact_cur.setZero();
    contact_duration.setZero();
    float phase_start_time(0.0);
    float phase_end_time(0.0);
    int phase_horizon(0.0);
    int n_phases(0);
    float t = 0.0;

    quad_ref_ptr->get_contact_at_t(contact_prev, t);
    quad_ref_ptr->get_contact_duration_at_t(contact_duration, t);

    while (approx_leq_scalar(t, plan_duration))
    {        
        // Get a contact at current time
        quad_ref_ptr->get_contact_at_t(contact_cur, t);       

        // Determine constructing one phase: change of contact or reach the planning horizon
        if ((contact_cur.cwiseNotEqual(contact_prev)).any() ||
            approx_geq_scalar(t, plan_duration))
        {                                
            phase_end_time = t;
            n_phases++;
            phase_horizon = (int)round((phase_end_time - phase_start_time) / dt_sim);

            pdata->phase_start_times.push_back(phase_start_time);
            pdata->phase_end_times.push_back(phase_end_time);
            pdata->phase_contacts.push_back(contact_prev);
            pdata->phase_horizons.push_back(phase_horizon);
            pdata->contact_durations.push_back(contact_duration);

            bool phase_reach_to_end = (contact_prev.cwiseNotEqual(contact_prev)).any();
            pdata->is_phase_reach_end.push_back(phase_reach_to_end);

            contact_prev = contact_cur;
            quad_ref_ptr->get_contact_duration_at_t(contact_duration, t);
            phase_start_time = phase_end_time;            
        }        
        pdata->n_phases = n_phases;
        t += dt_sim;
    }

    /* Initialize REB and AL parameters for ineq and eq constraints */
    const std::string& constraint_params_fname = "../HKDMPC/settings/constraint_params.info";
    loadConstrintParameters(constraint_params_fname, grf_reb_param, swing_reb_param,td_al_param);
    
    /* Create a phase and initialize its trajectory, cost, etc */
    for (int i = 0; i < n_phases; i++)
    {        
        shared_ptr<SinglePhase<T, 24, 24, 0>> phase;
        phase = make_shared<SinglePhase<T, 24, 24, 0>>();

        shared_ptr<Trajectory<T, 24, 24, 0>> traj;
        traj = make_shared<Trajectory<T, 24, 24, 0>>(dt_sim, pdata->phase_horizons[i]);

        /* Initialize state trajectory using the state reference */
        VecM<double, 24> xr_k;
        for (int k(0); k <= pdata->phase_horizons[i]; k++)
        {            
            hkd_reference.get_reference_at_t(xr_k, pdata->phase_start_times[i] + k * dt_sim);
            traj->X.at(k) = xr_k.cast<T>();
            traj->Xbar.at(k) = xr_k.cast<T>();
        }

        // Add trajectory to the phase (phase_idx)
        phase->set_trajectory(traj);

        create_problem_one_phase(phase, i);

        add_tconstr_one_phase(phase, i);

        phase->set_time_offset(pdata->phase_start_times[i] - pdata->phase_start_times[0]);

        phase->initialization();

        // Configure the set of shooting state
        phase->update_SS_config(pdata->phase_horizons[i] + 1);

        pdata->trajectory_ptrs.push_back(traj);
        pdata->phase_ptrs.push_back(phase);                
    }
    printf("\n");
    printf("Finished initializing HKDProbelm! \n\n");
}

/*
  @brief:   update at every dt_mpc
 */
template <typename T>
void HKDProblem<T>::update()
{
    /* update the multi-phase problem */
    for (int j = 0; j < nsteps_between_mpc; j++) // nsteps_between_mpc: number of simulation time steps between two mpc updates
    {
        // update the reference by one simulation time step
        quad_ref_ptr->step(dt_sim);

        float new_start_time = quad_ref_ptr->get_start_time();
        float new_end_time = quad_ref_ptr->get_end_time();

        // update the front end
        pdata->phase_start_times.front() += dt_sim;

        // check whether the first phase shrinks to a point
        if (approx_leq_scalar(pdata->phase_end_times.front(), new_start_time))
        {
            // If yes, remove the front phase
            pdata->pop_front_phase();
        }
        else
        {
            // else, pop_front one time step from the front phase and update phase horizon
            pdata->phase_ptrs.front()->pop_front();
            pdata->phase_horizons.front() --;
            pdata->phase_start_times.front() = new_start_time;
        }


        /* update the back end */
        
        VecM<int, 4> new_contact;
        quad_ref_ptr->get_contact_at_t(new_contact, new_end_time - new_start_time);
        bool contact_change = (new_contact.cwiseNotEqual(pdata->phase_contacts.back())).any();

        // If there is contact change, grow multi-phase problem by a new phase
        if (contact_change && pdata->is_phase_reach_end.back())
        {          
            float new_phase_start_time = pdata->phase_end_times.back();
            float new_phase_end_time = new_end_time;
            
            int new_phase_horizon = (int) round((new_phase_end_time - new_phase_start_time) / dt_sim);
            VecM<double, 4> new_contact_duration;
            quad_ref_ptr->get_contact_duration_at_t(new_contact_duration, new_end_time - new_start_time);
            pdata->phase_start_times.push_back(new_phase_start_time);
            pdata->phase_end_times.push_back(new_phase_end_time);            
            pdata->phase_horizons.push_back(new_phase_horizon);
            pdata->is_phase_reach_end.push_back(false);
            pdata->phase_contacts.push_back(new_contact);
            pdata->contact_durations.push_back(new_contact_duration);
            pdata->n_phases++;

            shared_ptr<Trajectory<T, 24, 24, 0>> traj_to_add;
            traj_to_add = make_shared<Trajectory<T, 24, 24, 0>>(dt_sim, new_phase_horizon);

            shared_ptr<SinglePhase<T, 24, 24, 0>> phase_to_add;
            phase_to_add = make_shared<SinglePhase<T, 24, 24, 0>>();

            phase_to_add->set_trajectory(traj_to_add);

            create_problem_one_phase(phase_to_add, pdata->n_phases - 1);

            phase_to_add->initialization();

            pdata->trajectory_ptrs.push_back(traj_to_add);

            pdata->phase_ptrs.push_back(phase_to_add);
        }
        // Else, grow the last phase by one time step
        else
        {
            pdata->phase_end_times.back() = new_end_time;
            pdata->phase_horizons.back()++;

            if (contact_change)
            {
                pdata->is_phase_reach_end.back() = true;
            }

            pdata->phase_ptrs.back()->push_back_default();
        }

        if (pdata->is_phase_reach_end.back())
        {
            add_tconstr_one_phase(pdata->phase_ptrs.back(), pdata->n_phases - 1);
        }
        
    }

    /* Update the multiple shooting configuration */
    for (int i = 0; i < pdata->n_phases; i++)
    {
        /* Update the time offset of each phase */
        pdata->phase_ptrs[i]->set_time_offset(pdata->phase_start_times[i] - pdata->phase_start_times[0]);

        pdata->phase_ptrs[i]->reset_params();

        if ((i == pdata->n_phases - 1 && pdata->phase_horizons[i] > 2) ||
            i < pdata->n_phases - 1)
        {
            pdata->phase_ptrs[i]->update_SS_config(pdata->phase_horizons[i] + 1);
        }

        pdata->trajectory_ptrs.front()->Ubar[0].setZero();
    }
}

template <typename T>
void HKDProblem<T>::create_problem_one_phase(shared_ptr<SinglePhase<T, 24, 24, 0>> phase, int idx)

{
    /* specialize dynamics and resetmap  */
    const auto &phase_contact = pdata->phase_contacts[idx];
    auto dynamics_callback = bind(&HKD::Model<T>::dynamics, &hkdModel,
                                  pc::_1, pc::_2, pc::_3, pc::_4, pc::_5, 
                                  phase_contact, (T)dt_sim);
    auto dynamics_partial_callback =
        bind(&HKD::Model<T>::dynamics_partial, &hkdModel,
             pc::_1, pc::_2, pc::_3, pc::_4, pc::_5,pc:: _6, pc::_7, 
             phase_contact, (T)dt_sim);

    /* set dynamics */
    phase->set_dynamics(dynamics_callback);
    phase->set_dynamics_partial(dynamics_partial_callback);

    /* set tracking cost */
    shared_ptr<HKDTrackingCost<T>> track_cost;
    track_cost = make_shared<HKDTrackingCost<T>>(phase_contact);
    track_cost->set_reference(&hkd_reference);
    phase->add_cost(track_cost);

    /* Set foot regularization */
    shared_ptr<HKDFootPlaceReg<T>> foot_reg;
    foot_reg = make_shared<HKDFootPlaceReg<T>>(phase_contact);
    foot_reg->set_quad_reference(quad_ref_ptr);
    phase->add_cost(foot_reg);

    /* Set GRF constraints if any*/
    if (phase_contact.cwiseEqual(1).any())
    {
        shared_ptr<GRFConstraint<T>> grfConstraint;
        grfConstraint = std::make_shared<GRFConstraint<T>>(phase_contact);
        grfConstraint->update_horizon_len(pdata->phase_horizons[idx]);
        grfConstraint->create_data();
        grfConstraint->initialize_params(grf_reb_param);
        phase->add_pathConstraint(grfConstraint);
    }

}

template <typename T>
void HKDProblem<T>::add_tconstr_one_phase(shared_ptr<SinglePhase<T, 24, 24, 0>> phase, int idx)
{
    /* Determine the touchdown status */
    const VecM<int, 4> &phase_contact_cur = pdata->phase_contacts[idx];
    VecM<int, 4> touchdown_status;
    VecM<int, 4> phase_contact_next;
    touchdown_status.setZero();
    if (idx < pdata->n_phases - 1) // If it is an intermediate phase
    {
        phase_contact_next = pdata->phase_contacts[idx + 1];
    }
    else // else the last phase
    {
        quad_ref_ptr->get_contact_at_t(phase_contact_next, plan_duration + dt_mpc);
    }

    for (int leg = 0; leg < 4; leg++)
    {
        if (phase_contact_cur[leg] == 0 && phase_contact_next[leg] == 1)
        {
            touchdown_status[leg] = 1;
        }
    }

    /* set resetmap */
    auto resetmap_callback = bind(&HKDReset<T>::resetmap, &hkdReset,
                                  pc::_1, pc::_2, phase_contact_cur, phase_contact_next);
    auto resetmap_partial_callback = bind(&HKDReset<T>::resetmap_partial, &hkdReset,
                                          pc::_1, pc::_2, phase_contact_cur, phase_contact_next);

    phase->set_resetmap(resetmap_callback);
    phase->set_resetmap_partial(resetmap_partial_callback);

    /* Set touchdown constraint if any */
    if (find_eigen(touchdown_status, 1).size() > 0)
    {
        shared_ptr<TouchDownConstraint<T>> tdConstraint;
        tdConstraint = std::make_shared<TouchDownConstraint<T>>(touchdown_status);
        tdConstraint->create_data();
        tdConstraint->initialize_params(td_al_param);
        phase->add_terminalConstraint(tdConstraint);
    }
}


template <typename T>
void HKDProblem<T>::pretty_print()
{
    tabulate::Table config_print;
    config_print.add_row({"Plan duration", "dt_sim", "dt_mpc"});
    config_print.add_row({std::to_string(plan_duration),
                          std::to_string(dt_sim),
                          std::to_string(dt_mpc)});

    config_print.column(1).format() 
        .font_align(tabulate::FontAlign::center);

    for (size_t i = 0; i < config_print[0].size(); i++)
    {
        config_print[0][i].format()
        .font_color(tabulate::Color::yellow)
        .font_align(tabulate::FontAlign::center)
        .font_style({tabulate::FontStyle::bold});
    }

    std::cout << std::endl << config_print << std::endl;
                             
    tabulate::Table problem;
    problem.add_row({"Phase Index", "Horizon", "Start Time", "End Time", "Contact", "Contact Duration"});    
    for (int i = 0; i < pdata->n_phases; i++)
    {
        problem.add_row({std::to_string(i), std::to_string(pdata->phase_horizons[i]),
                         std::to_string(pdata->phase_start_times[i]), 
                         std::to_string(pdata->phase_end_times[i]),
                         eigenToString(pdata->phase_contacts[i].transpose()),
                         eigenToString(pdata->contact_durations[i].transpose())});                            
    }    

    // center-align and color header cells
    for (size_t i = 0; i < problem[0].size(); ++i) {
        problem[0][i].format()
        .font_color(tabulate::Color::yellow)
        .font_align(tabulate::FontAlign::center)
        .font_style({tabulate::FontStyle::bold});

         problem.column(i).format()
            .font_align(tabulate::FontAlign::center);       
    }    

    std::cout << problem << std::endl << std::endl;
}

template <typename T>
void HKDProblem<T>::lcm_publish()
{
    // int k(0);
    // float timeStart = ref_data->startTimes[0];

    // reset_lcm_data();

    // for (int i(0); i < ref_data->n_phases; i++)
    // {
    //     for (int j(0); j < ref_data->horizons[i]; j++)
    //     {
    //         lcm_pdata.times.push_back((timeStart + k * (ref_data->dt)) * 1e06);

    //         for (int l(0); l < 4; l++)
    //         {
    //             lcm_pdata.contacts[l].push_back(ref_data->contactSeq[i][l]);
    //             lcm_pdata.qdummy_r[3 * l].push_back(ref_data->Xr[i][j][12 + 3 * l]);
    //             lcm_pdata.qdummy_r[3 * l + 1].push_back(ref_data->Xr[i][j][12 + 3 * l + 1]);
    //             lcm_pdata.qdummy_r[3 * l + 2].push_back(ref_data->Xr[i][j][12 + 3 * l + 2]);

    //             lcm_pdata.qdummy[3 * l].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3 * l]);
    //             lcm_pdata.qdummy[3 * l + 1].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3 * l + 1]);
    //             lcm_pdata.qdummy[3 * l + 2].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3 * l + 2]);
    //         }

    //         for (int d = 0; d < 3; d++)
    //         {
    //             lcm_pdata.eul_r[d].push_back(ref_data->Xr[i][j][d]);
    //             lcm_pdata.pos_r[d].push_back(ref_data->Xr[i][j][3 + d]);
    //             lcm_pdata.vel_r[d].push_back(ref_data->Xr[i][j][6 + d]);
    //             lcm_pdata.omega_r[d].push_back(ref_data->Xr[i][j][9 + d]);

    //             lcm_pdata.eul[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][d]);
    //             lcm_pdata.pos[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][3 + d]);
    //             lcm_pdata.vel[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][6 + d]);
    //             lcm_pdata.omega[d].push_back(pdata->trajectory_ptrs[i]->Xbar[j][9 + d]);
    //         }
    //         k++;
    //     }
    // }

    // lcm_pdata.n_timesteps = k;
    // _lcm.publish("DEBUG_HKDMPC", &lcm_pdata);
}

template <typename T>
void HKDProblem<T>::reset_lcm_data()
{
    // lcm_pdata.times.clear();
    // lcm_pdata.contacts.clear();
    // lcm_pdata.pos_r.clear();
    // lcm_pdata.eul_r.clear();
    // lcm_pdata.vel_r.clear();
    // lcm_pdata.omega_r.clear();
    // lcm_pdata.qdummy_r.clear();

    // lcm_pdata.pos.clear();
    // lcm_pdata.eul.clear();
    // lcm_pdata.vel.clear();
    // lcm_pdata.omega.clear();
    // lcm_pdata.qdummy.clear();

    // lcm_pdata.contacts.resize(4);
    // lcm_pdata.pos_r.resize(3);
    // lcm_pdata.eul_r.resize(3);
    // lcm_pdata.vel_r.resize(3);
    // lcm_pdata.omega_r.resize(3);
    // lcm_pdata.qdummy_r.resize(12);
    // lcm_pdata.pos.resize(3);
    // lcm_pdata.eul.resize(3);
    // lcm_pdata.vel.resize(3);
    // lcm_pdata.omega.resize(3);
    // lcm_pdata.qdummy.resize(12);
}

template class HKDProblem<double>;