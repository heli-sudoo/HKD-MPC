#include "HSDDP_Utils.h"

std::ostream &print_time(std::ostream &os, std::vector<TIME_PER_ITERATION> &time_ddp)
{
    for (size_t i = 0; i < time_ddp.size(); i++)
    {
        os << time_ddp[i].DDP_iter << "\t"
           << time_ddp[i].n_bws << "\t" << time_ddp[i].time_bws << "\t"
           << time_ddp[i].n_fit << "\t" << time_ddp[i].time_fit << "\t"
           << time_ddp[i].time_partial << std::endl;
    }
    return os;
}