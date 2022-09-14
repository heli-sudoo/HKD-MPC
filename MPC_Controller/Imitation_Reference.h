#ifndef IMITATION_REFERENCE_H
#define IMITATION_REFERENCE_H

#include <deque>
#include <string>
#include "HSDDP_CPPTypes.h"
#include "HKDReference.h"

using std::deque;
using std::string;

template<typename T>
class Imitation_Reference
{
private:    
    HKDReferenceData<T> data;

public:
    Imitation_Reference(){}

    void load_state_data(const string&);

    void load_contact_data(const string&);

    void compute_status_duration();

    HKDReferenceData<T>* get_data_ptr() {return &data;}
};



#endif