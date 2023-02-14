#ifndef HSDDP_COMPOUNDTYPES_H
#define HSDDP_COMPOUNDTYPES_H

#include "HSDDP_CPPTypes.h"
#include <cstring>

// Structure collecting all HSDDP parameters
struct HSDDP_OPTION
{
    double alpha = 0.1;               // line search udpate paramer
    double gamma = 0.01;              // scale the expected cost reduction
    double update_penalty = 8;        // penalty update parameter
    double update_relax = 0.1;        // relaxation parameter udpate parameter
    double update_regularization = 2; // regularization parameter update parameter
    double update_ReB = 7;            // update barrier function weighting
    double max_DDP_iter = 3;          // maximum inner loop iteration
    double max_AL_iter = 2;           // maximum outer loop iteration/*  */
    double DDP_thresh = 1e-03;        // inner loop convergence threshhold
    double tconstr_thresh = 1e-03;
    double pconstr_thresh = 1e-03;
    bool AL_active = 1;               // activate terminal constraint
    bool ReB_active = 1;              // activate path constraint
    bool smooth_active = 0;           // activate control smoothness penalization
    std::vector<int> SS_set;          // shooting state set
};


#ifdef TIME_BENCHMARK
extern vector<TIME_PER_ITERATION> time_ddp;
#endif //TIME_BENCHMARK

#endif // HSDDP_COMPOUNDTYPES_H