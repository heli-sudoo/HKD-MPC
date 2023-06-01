/*!
 * @file opt_main.cpp
 * @brief Main Function for the standalone DDP trajectory optimizer
 *
 * The main function initilizes the DDP solver updates upon the request
 * of the main Imitation Controller
 */

#include "HKDMPC.h"
#include "HSDDP_Utils.h"

#ifdef TIME_BENCHMARK
#include <chrono>
using namespace std::chrono;
using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;
static int fit_iter = 0;
#endif // TIME_BENCHMARK

template <typename T>
void HKDMPCSolver<T>::initialize()
{
    printf("Initializing HKDMPC Solver ... \n\n");
    std::string fname_ddp_setting("../HKDMPC/settings/ddp_setting.info");
    loadHSDDPSetting(fname_ddp_setting, ddp_options);

    mpc_config.plan_duration = .6;
    mpc_config.nsteps_between_mpc = 1;
    mpc_config.timeStep = 0.01;
    dt_mpc = mpc_config.timeStep;

    opt_problem.clear_problem_data();
    opt_problem_data.quad_ref_ptr = &quad_reference;
    opt_problem.set_problem_data(&opt_problem_data, mpc_config);
    opt_problem.initialization();

#ifdef DEBUG_MODE
    opt_problem.pretty_print();
#endif    

    mpc_time = 0;
    mpc_time_prev = 0;

    // set the initial condition
    xinit.setZero(24);
    body << 0, 0, 0, 0, 0, 0.2486, 0, 0, 0, 0, 0, 0;
    qJ << 0, -0.8, 1.6, 0, -0.8, 1.6, 0, -0.8, 1.6, 0, -0.8, 1.6;
    pos = body.segment(3, 3);
    eul = body.head(3);

    // get the contact status of the very first phase
    const auto &initial_contact = opt_problem_data.phase_contacts.front();
    compute_hkd_state(eul, pos, qJ, qdummy, initial_contact);

    xinit << body, qdummy;

    // Assemble the multi-phase probelm and solve it with MSDDP
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (auto phase : opt_problem_data.phase_ptrs)
    {
        multiple_phases.push_back(phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(xinit);

#ifdef TIME_BENCHMARK
    auto start = high_resolution_clock::now();
#endif
    solver.solve(ddp_options);
#ifdef TIME_BENCHMARK
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
    solve_time = duration.count();
#endif
#ifdef DEBUG_MODE
    solver.get_solver_info(solver_info.cost,
                           solver_info.dyn_feas,
                           solver_info.eqn_feas,
                           solver_info.ineq_feas);
    publish_solver_info();    
#endif    
    std::string folder_name = "../HKDMPC/log/";
    log_trajectory_sequence(folder_name, opt_problem_data.trajectory_ptrs);
    mpc_iter = 0;

    printf("Finished initializing HKDMPC Solver \n\n");

    update_foot_placement();
    publish_mpc_cmd();
    
#ifdef DEBUG_MODE    
    publish_debugfoot();    
#endif
}

template <typename T>
void HKDMPCSolver<T>::update()
{
    mpc_mutex.lock(); // lock mpc to prevent updating while the previous hasn't finished

    // use less iterations when resolving DDP
    ddp_options.max_AL_iter = 2;
    ddp_options.max_DDP_iter = 1;
    mpc_iter++;

    printf("************************************* \n");
    printf("************Resolving MPC************ \n");
    printf("********MPC Iteration = %d*********** \n", mpc_iter);
    printf("************************************* \n");

    /* update the problem */
    opt_problem.update();

#ifdef DEBUG_MODE    
    opt_problem.pretty_print();
#endif

    /* update current state*/
    eul << hkd_data.rpy[2], hkd_data.rpy[1], hkd_data.rpy[0];
    pos << hkd_data.p[0], hkd_data.p[1], hkd_data.p[2];
    vel << hkd_data.vWorld[0], hkd_data.vWorld[1], hkd_data.vWorld[2];
    omega << hkd_data.omegaBody[0], hkd_data.omegaBody[1], hkd_data.omegaBody[2];
    for (int i(0); i < 12; i++)
    {
        qJ[i] = hkd_data.qJ[i];
    }
    const auto &current_contact = opt_problem_data.phase_contacts.front();
    compute_hkd_state(eul, pos, qJ, qdummy, current_contact);
    xinit << eul, pos, omega, vel, qdummy;

    /* build solver and solve the TO problem */
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (auto phase : opt_problem_data.phase_ptrs)
    {
        multiple_phases.push_back(phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(xinit);
#ifdef TIME_BENCHMARK
    auto start = high_resolution_clock::now();
#endif
    solver.solve(ddp_options);
#ifdef TIME_BENCHMARK
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
    solve_time = duration.count();
    printf("solve time = %f \n", solve_time);
#endif

#ifdef DEBUG_MODE
    solver.get_solver_info(solver_info.cost,
                           solver_info.dyn_feas,
                           solver_info.eqn_feas,
                           solver_info.ineq_feas);
    publish_solver_info();
#endif

    update_foot_placement();
    publish_mpc_cmd();
#ifdef DEBUG_MODE
    // publish_debugfoot();
    // opt_problem.lcm_publish();
#endif
    mpc_mutex.unlock();
}

template <typename T>
void HKDMPCSolver<T>::mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                                          const hkd_data_lcmt *msg)
{
    (void) (rbuf);
    (void) (chan);
    mpc_mutex.lock();
    printf(GRN);
    printf("Received resolving request\n");
    printf(RESET);

    if (msg->reset_mpc)
    {
        ddp_options.MS = msg->MS;
        mpc_mutex.unlock();
        initialize();
        return;
    }

    std::memcpy(&hkd_data, msg, sizeof(hkd_data));
    mpc_time_prev = mpc_time;
    mpc_time = hkd_data.mpctime;

    // get the current foot placements
    const auto &current_pf = hkd_data.foot_placements;
    for (int l = 0; l < 4; l++)
    {
        pf[l] << current_pf[3 * l], current_pf[3 * l + 1], current_pf[3 * l + 2];
    }
    mpc_mutex.unlock();
    std::thread solve_mpc_thread(&HKDMPCSolver::update, this);
    solve_mpc_thread.detach(); // detach the thread from the main thread. The thread would exit once it completes
}

/*
    @brief: Go through the trajectory to find the next foot placement
            Needs to update every time after an MPC update
*/
template <typename T>
void HKDMPCSolver<T>::update_foot_placement()
{
    // mpc_mutex.lock();
    int found_next_ctact[4] = {0};
    const int &n_phases = opt_problem_data.n_phases;
    const auto &ctactSeq = opt_problem_data.phase_contacts;

    for (int i(0); i < n_phases - 1; i++)
    {
        const auto &ctact = ctactSeq[i];
        const auto &ctactn = ctactSeq[i + 1];
        for (int l(0); l < 4; l++)
        {
            // If we havn't found the next contact
            if (!found_next_ctact[l])
            {
                // Check the search pattern [0, 1]. Iterate through the contact sequence until match the search
                // pattern. qdummy associated with 1 and leg l is then foot placement
                if (ctact[l] == 0 && ctactn[l] == 1)
                {
                    const auto& qdummy_ = opt_problem_data.trajectory_ptrs[i + 1]->Xbar[0].tail(12);
                    pf[l] = qdummy_.segment(3 * l, 3).template cast<float>();
                    found_next_ctact[l] = 1;
                }
            }
        }
        // Break if not found after four phases
        if (i >= 4)
        {
            break;
        }
    }
    // mpc_mutex.unlock();
}

template <typename T>
void HKDMPCSolver<T>::publish_mpc_cmd()
{
    int num_controls = mpc_config.nsteps_between_mpc;
    num_controls += 7; // use 7 more controls than control duration to account for delay

    hkd_cmds.N_mpcsteps = num_controls;
    auto &trajseq = opt_problem_data.trajectory_ptrs;
    auto &ctactSeq = opt_problem_data.phase_contacts;
    auto &statusDuration = opt_problem_data.contact_durations;
    int k(0), s(0), i(0);

    while (k < hkd_cmds.N_mpcsteps)
    {
        if (s >= trajseq[i]->horizon)
        {
            s = 0;
            i++;
        }
        for (int j = 0; j < 24; j++)
        {
            hkd_cmds.hkd_controls[k][j] = trajseq[i]->Ubar[s][j];
        }
        for (int j = 0; j < 12; j++)
        {
            hkd_cmds.des_body_state[k][j] = trajseq[i]->Xbar[s][j];
        }
        for (int m = 0; m < 12; m++)
        {
            for (int n = 0; n < 12; n++)
            {
                hkd_cmds.feedback[k][m][n] = trajseq[i]->K[s](m, n);
            }
        }

        hkd_cmds.mpc_times[k] = mpc_time + (k * dt_mpc);
        for (int l = 0; l < 4; l++)
        {
            hkd_cmds.contacts[k][l] = ctactSeq[i][l];
            hkd_cmds.statusTimes[k][l] = statusDuration[i][l];
        }

        s++;
        k++;
    }
    for (int l = 0; l < 4; l++)
    {
        hkd_cmds.foot_placement[3 * l] = pf[l][0];
        hkd_cmds.foot_placement[3 * l + 1] = pf[l][1];
        hkd_cmds.foot_placement[3 * l + 2] = pf[l][2];
    }
    hkd_cmds.solve_time = solve_time;
    mpc_lcm.publish("mpc_command", &hkd_cmds);
    printf(GRN);
    printf("published a mpc command message \n");
    printf(RESET);
}

template <typename T>
void HKDMPCSolver<T>::publish_debugfoot()
{
    debug_foot_data.N = 0;
    debug_foot_data.contacts.clear();
    debug_foot_data.qdummy.clear();
    const auto &contacts = opt_problem_data.phase_contacts;
    int n_phases = opt_problem_data.n_phases;
    for (int i = 0; i < n_phases; i++)
    {
        auto traj = opt_problem_data.trajectory_ptrs[i];
        auto horizon = traj->horizon;
        for (int k = 0; k < horizon; k++)
        {
            vector<float> qdummy_(traj->Xbar[k].data() + 12, traj->Xbar[k].data() + 23);
            vector<int> ctact(contacts[i].data(), contacts[i].data() + 3);
            debug_foot_data.contacts.push_back(ctact);
            debug_foot_data.qdummy.push_back(qdummy_);
        }
        debug_foot_data.N += horizon;
    }
    mpc_lcm.publish("debug_foot", &debug_foot_data);
}

template <typename T>
void HKDMPCSolver<T>::publish_solver_info()
{
    solver_info.n_iter = solver_info.cost.size();
    solver_info_lcm.publish("solver_info", &solver_info);
}

template class HKDMPCSolver<double>;