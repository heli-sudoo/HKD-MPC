#ifndef HSDDP_UTILS_H
#define HSDDP_UTILS_H
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


using std::vector;
using std::ostream;

struct  TIME_PER_ITERATION
{
    int DDP_iter = 0;
    int n_bws = 0; // number of backward sweep iterations
    double time_bws =0; // backward sweep time
    int n_fit = 0; // number of iterations in line search
    double time_fit = 0; // line search time
    double time_partial = 0; // time for computing derivatives (dynamics and cost)

};

template<typename T>
void unused_ignore(T &&){} // && reference to rvalue (C++ 11)

template<typename Derived, typename T>
vector<int>find_eigen(const Eigen::DenseBase<Derived>& v, const T& val){
    typedef typename Derived::Scalar val_type;
    vector<int> index;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == static_cast<val_type>(val))
        {
            index.push_back(i);
        }
        
    }
    return index;
}

template<typename SeqContainer>
void append_vector(SeqContainer& v1, const SeqContainer& v2){
    for (auto &val : v2)
    {
        v1.push_back(val);
    }    
}

// Check if two floating-point number are almost equal
template<typename T>
bool almostEqual_number(T n1, T n2)
{
    T tol = 1e-8;
    T err = std::abs(n1 - n2);
    if (err <= tol)
    {
        return true;
    }
    return false;
}

ostream &print_time(ostream &os, std::vector<TIME_PER_ITERATION> &time_ddp);


#endif