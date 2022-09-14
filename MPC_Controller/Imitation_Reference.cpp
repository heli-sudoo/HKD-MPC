#include <fstream>
#include <iostream>
#include "Imitation_Reference.h"
#include "cTypes.h" // colored print

using std::fstream;
using std::stringstream;

template <typename T>
void Imitation_Reference<T>::load_contact_data(const string &ctact_fname)
{
    // clear contact state information
    data.clear();
    int &n_phases = data.n_phases;

    // read contact state information from file
    fstream fstrm(ctact_fname);
    string line;
    if (!fstrm.is_open())
    {
        printf("Failed to open contact file \n");
        return;
    }
    printf("Reading reference contact information \n");

    // Get the first line and print the names of each column
    getline(fstrm, line);
    std::cout << line << std::endl;

    // Go through the rest of lines
    T sT, eT;
    int horizon;
    VecM<int, 4> contact;
    string word;
    n_phases = 0;
    while (getline(fstrm, line))
    {
        stringstream s(line); // break line to s with delimeter (,)
        // std::cout << line << std::endl;
        int i = 0;
        while (getline(s, word, ','))
        {
            if (i <= 3)
            {
                contact[i] = std::stoi(word);
            }
            if (i == 4)
            {
                sT = static_cast<T>(std::stod(word));
            }
            if (i == 5)
            {
                eT = static_cast<T>(std::stod(word));
            }
            if (i == 6)
            {
                horizon = std::stoi(word);
            }

            i++;
        }
        data.contactSeq.push_back(contact);
        data.startTimes.push_back(sT);
        data.endTimes.push_back(eT);
        data.horizons.push_back(horizon);
        n_phases++;
    }
    fstrm.close();
    long int eT_ns = 1e9 * eT;
    long int sT_ns = 1e9 * sT;
    data.dt = ((T)(eT_ns - sT_ns)) / (1e9 * horizon);
}

template <typename T>
void Imitation_Reference<T>::load_state_data(const string &state_fname)
{
    auto &Xr = data.Xr;
    auto &Ur = data.Ur;
    auto &Yr = data.Yr;
    auto &horizons = data.horizons;

    // read reference state trajectory data from file
    fstream fstrm(state_fname);
    if (!fstrm.is_open())
    {
        printf("Failed to open state file \n");
        return;
    }
    printf("Reading reference state data \n");

    VecM<T, 24> x;
    VecM<T, 24> u;
    VecM<T, 0> y;
    x.setZero();
    u.setZero();
    y.setZero();

    deque<VecM<T, 24>> Xr_phase;
    deque<VecM<T, 24>> Ur_phase;
    deque<VecM<T, 0>> Yr_phase;

    string line, word;
    int count(0);
    int pidx(0);
    while (getline(fstrm, line))
    {
        stringstream s(line);
        int i(0), idx(0);
        while (getline(s, word, ','))
        {
            (i == 0) ? idx = std::stoi(word) : x[i - 1] = static_cast<T>(std::stod(word));
            i++;
        }
        if (idx == pidx)
        {
            Xr_phase.push_back(x);
        }
        else // idx > pidx
        {
            Xr.push_back(Xr_phase);
            Xr_phase.clear();
            Xr_phase.push_back(x);
            pidx++;
        }
    }
    fstrm.close();
    printf("finished reading %i phases of state trajectory \n", pidx);

    // check size
    if (Xr.size() != data.contactSeq.size() - 1)
    {
        printf(RED);
        printf("number of contact phases does not match contact file \n");
        printf(RESET);
        return;
    }

    // assemble reference control and output
    for (int i = 0; i < Xr.size(); i++)
    {
        // check size
        if (Xr[i].size() != horizons[i] + 1)
        {
            printf(RED);
            printf("state trajectory size does not match contact file at phase %i\n", i);
            printf("trajectory size = %lu \n", Xr[i].size());
            printf("horizon = %i \n", horizons[i]);
            printf(RESET);
            return;
        }
        u.setZero();
        for (int l(0); l < 4; l++)
        {
            if (data.contactSeq[i][l])
            {
                u[3*l+2] = 90/data.contactSeq[i].sum();
            }
            
        }
        
        Ur_phase.clear();
        Yr_phase.clear();
        Ur_phase.assign(horizons[i], u);
        Yr_phase.assign(horizons[i], y);
        Ur.push_back(Ur_phase);
        Yr.push_back(Yr_phase);
    }
    printf("finished loading imimitation reference \n \n");
}

template <typename T>
void Imitation_Reference<T>::compute_status_duration()
{
    data.statusDuration.setZero(4, data.contactSeq.size());
    for (int l(0); l < 4; l++)
    {
        int idx_s(0), idx_e(0), count(0);
        int ctact_s = data.contactSeq[0][l];
        for (int i(1); i < data.contactSeq.size(); i++)
        {
            const int &ctact_i = data.contactSeq[i][l];
            // if contact status does not changed, increment idx_e
            if (ctact_i == ctact_s)
            {
                idx_e++;
            }
            else
            {
                count = idx_e - idx_s + 1;
                T dur = data.endTimes[idx_e] - data.startTimes[idx_s];
                data.statusDuration.row(l).middleCols(idx_s, count) = DVec<T>::Constant(count, dur).transpose();

                ctact_s = ctact_i;
                idx_s = i;
                idx_e = i;
            }
        }
        // do some wrapup for the last phase
        count = idx_e - idx_s + 1;
        T dur = data.endTimes[idx_e] - data.startTimes[idx_s];
        data.statusDuration.row(l).middleCols(idx_s, count) = DVec<T>::Constant(count, dur).transpose();
    }
}

template class Imitation_Reference<double>;