#include <fstream>
#include <iostream>
#include "QuadReference.h"
#include "HSDDP_Utils.h"

void QuadReference::initialize(float plan_horizon)
{
    printf("Initializing reference trajectory ... \n");

    data.clear();
    t_cur = 0;
    k_cur = 0;
    dt = tp_data.dt;
    dur = plan_horizon;    
    data.start_time = t_cur;
    data.end_time = t_cur + dur;
    sz = (int)round(plan_horizon / dt) + 1;

    std::cout << "copying " << sz + 1 << " data (unit: time step)" << std::endl;

    // Copy one more elements greater than the planning horizon
    // Used to reason about resetmap and terminal constraints
    std::copy(tp_data.begin(), tp_data.begin() + sz + 1, std::back_inserter(data.get_container())); 

    printf("Finished initializing reference trajectory! \n\n");
}

/*
    @brief: Shift the reference data forward by one time step (simulation time step)
    @params: for hkd mpc, we use dt_sim = 10 ms. Currently only suports dt_ref <= dt_sim
    To Do:  consider the case where dt_sim < dt_ref (whole-body dynamics)
*/
void QuadReference::step(float dt_sim)
{
    for (int i = 1; approx_leq_scalar(i * dt, dt_sim); i++)
    {
        k_cur++;
        t_cur += dt;

        data.pop_front();        
        data.push_back(tp_data[k_cur + sz]);

        // shifting the reference time
        data.start_time = t_cur;
        data.end_time = t_cur + dur;
    }
}

/* 
    @brief: Update the foot placements based on the current body state using Raibert Heuristics
*/
void QuadReference::update_foot_placements(const VecM<double, 12>& body_state_t)
{
    // To Be Implemented
    (void) (body_state_t);
}

/*
    brief: query a reference at time t
    params:
        t:  Relative time with repect to the current reference (always starts at zero)
            Not necessarily a multiple of dt. Can be any time between/on two adjcent states
    return: Pointer to the reference state at time t
*/
QuadAugmentedState *QuadReference::get_a_reference_ptr_at_t(float t)
{
    int k = (int)floor(t / dt); // clip the time index towards left

    if (t - k*dt > 0.5*dt) k ++; //avoid floating point precision problem
    

    // If queried time reaches the end, use the last element
    if (k > sz)
    {
        printf("warning: queried reference out of scope \n");
        printf("queried time = %f \n\n", t + t_cur);
        k = sz;
    }
    return data.get_ptr(k);
}


/*
    @brief: Get contact status at time t
*/
void QuadReference::get_contact_at_t(VecM<int, 4> &contact, float t)
{    
    int k = (int)floor(t / dt);

    if (t - k*dt > 0.5*dt) k ++; //avoid floating point precision problem

    if (k > sz)
    {        
        printf("warning: queried reference out of scope \n");
        printf("queried time = %f \n\n", t + t_cur);        
        k = sz;
    }
    
    contact = data[k].contact;
}

/*
    @brief: Get the duration of contact status at time t
*/
void QuadReference::get_contact_duration_at_t(VecM<double, 4>& contact_dur, float t)
{
    int k = (int)floor(t / dt);

    if (t - k*dt > 0.5*dt) k ++; //avoid floating point precision problem

    if (k > sz)
    {
        printf("warning: queried reference out of scope \n");
        printf("queried time at %f \n\n", t + t_cur);
        k = sz;
    }
    
    contact_dur = data[k].status_dur;
}

/*
    @brief: Load reference trajectory from txt file
    @parameters:
            fname: path/name of the reference trajectory file
            reorder: re-order the left and right legs
                     HKD-MPC uses the same convention as in Cheetah Software
                     MHPC employs the convention as in the urdf file
*/
void QuadReference::load_top_level_data(const std::string& fname, bool reorder)
{
    printf("Loading reference trajectory file ...\n");    

    std::fstream fstrm(fname, std::fstream::in);
    if (!fstrm.is_open())
    {
        std::cout <<"Error opening the file " << fname << std::endl;
        return;
    }
    std::cout << "Openning file " << fname << std::endl;

    int bs_sz(12), nLegs(4);
    std::string line, word;
    QuadAugmentedState quad_state;
    quad_state.SetZero();
    

    // Loop through each line of the file
    while (getline(fstrm, line))
    {        
        if (line == "dt")
        {            
            getline(fstrm, line);           // Get the next line
            tp_data.dt = std::stof(line);

            continue;
        }
        
        if (line.find("body_state") != std::string::npos)
        {
            quad_state.SetZero();
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.body_state[i] = std::stof(word);
                i++;
                if (i >= bs_sz) break;                
            }            
            
            continue;  
        }    

        if (line.find("qJ") != std::string::npos)
        {
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {                
                quad_state.qJ[i] = std::stof(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            continue;
        }       
       
        if (line.find("foot_placements") != std::string::npos)
        {
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.foot_placements[i] = std::stof(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            continue;
        }       

        if (line.find("grf") != std::string::npos)
        {
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.grf[i] = std::stof(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            continue;
        }        

        if (line.find("torque") != std::string::npos)
        {
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.torque[i] = std::stof(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            continue;
        }       

        if (line.find("contact") != std::string::npos)
        {
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.contact[i] = std::stoi(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            continue;
        }       

        if (line.find("status_dur") != std::string::npos)
        {            
            getline(fstrm, line);           // Get the next line
            std::stringstream lstrm(line);  // Break the line to words (default delimeter is " ")

            int i = 0;
            while (lstrm >> word)
            {
                quad_state.status_dur[i] = std::stof(word);
                i++;
                if (i >= 3 * nLegs) break;
            }    

            tp_data.push_back(quad_state);           
        }       
    }
    
    fstrm.close();

    /* Print the information of loaded trajectory */
    printf("Finished loading data \n");
    printf("Loaded trajectory size = %ld \n\n", tp_data.size());

    /* Reorder the reference state if needed (used for MHPC)*/
    if (reorder)
    {
        reorder_states();

    }
   
}

void QuadReference::reorder_states()
{
    QuadAugmentedState state_reordered;
    for (size_t i(0); i < tp_data.size(); ++i)
    {   
        const QuadAugmentedState& state_org = tp_data[i];

        /* reorder body state. original: [eul, pos, omega, vWorld]. reordered: [pos, eul, vWorld, eul_rate] */
        state_reordered.body_state << state_org.body_state.segment<3>(3), state_org.body_state.head<3>(),
                                      state_org.body_state.tail<3>(), state_org.body_state.segment<3>(6);
        state_reordered.body_state[2] = 0.25;       
        // state_reordered.body_state[6] = 1.0;                                      
        /* flip the right and left legs */
        state_reordered.qJ << state_org.qJ.segment<3>(3), state_org.qJ.head<3>(), state_org.qJ.tail<3>(), state_org.qJ.segment<3>(6);       
        // state_reordered.qJ << Vec3<double>(0, -0.8, 1.6).replicate<4,1>();

        state_reordered.qJd.setZero();
        state_reordered.foot_placements << state_org.foot_placements.segment<3>(3), state_org.foot_placements.head<3>(), 
                                           state_org.foot_placements.tail<3>(), state_org.foot_placements.segment<3>(6);

        state_reordered.grf << state_org.grf.segment<3>(3), state_org.grf.head<3>(), state_org.grf.tail<3>(), state_org.grf.segment<3>(6);                                               

        state_reordered.torque << state_org.torque.segment<3>(3), state_org.torque.head<3>(), state_org.torque.tail<3>(), state_org.torque.segment<3>(6);

        state_reordered.contact << state_org.contact[1], state_org.contact[0], state_org.contact[3], state_org.contact[2];        
        // state_reordered.contact.setOnes();
        state_reordered.status_dur << state_org.status_dur[1], state_org.status_dur[0], state_org.status_dur[3], state_org.status_dur[2];

        state_reordered.qJ(Eigen::seqN(1,4,3)) = -state_reordered.qJ(Eigen::seqN(1,4,3)); 
        state_reordered.qJ(Eigen::seqN(2,4,3)) = -state_reordered.qJ(Eigen::seqN(2,4,3));         
        state_reordered.torque(Eigen::seqN(1,4,3)) = -state_reordered.torque(Eigen::seqN(1,4,3)); 
        state_reordered.torque(Eigen::seqN(2,4,3)) = -state_reordered.torque(Eigen::seqN(2,4,3)); 

        tp_data[i] = state_reordered;
    }
    
}


