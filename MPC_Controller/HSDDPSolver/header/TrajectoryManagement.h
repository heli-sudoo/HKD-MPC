#ifndef TRAJECTORY_MANAGEMENT_H
#define TRAJECTORY_MANAGEMENT_H

#include <deque>
#include <typeinfo>
#include <memory> // shared_ptr (since C++ 11)
#include <vector>
#include <deque>
#include <type_traits> // hearder for is_scalar()
#include "HSDDP_CPPTypes.h"
#include "CostBase.h"

using std::vector;
using std::shared_ptr;
using std::deque;

using namespace std;

template <class Model_>
struct BriefTrajectoryIterator
{
    typedef typename Model_::Scalar T;
    static const int xs = Model_::state_dim;
    static const int us = Model_::contrl_dim;
    static const int ys = Model_::output_dim;

    typename deque<VecM<T, xs>>::iterator Xiter;
    typename deque<VecM<T, us>>::iterator Uiter;
    typename deque<VecM<T, ys>>::iterator Yiter;
};

template <class Model_>
struct BriefTrajectory
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef typename Model_::Scalar T;
    static const int xs = Model_::state_dim;
    static const int us = Model_::contrl_dim;
    static const int ys = Model_::output_dim;

    deque<VecM<T, xs>> X;
    deque<VecM<T, us>> U;
    deque<VecM<T, ys>> Y;
    size_t len;
    size_t length()
    {
        // check_length();
        len = X.size();
        return len;
    }
    void check_length()
    {
        if (X.size() != U.size())
        {
            printf("Trajectory length of state and control do not match \n");
        }
        if (X.size() != Y.size())
        {
            printf("Trajectory length of state and output do not match \n");
        }
    }
    void clear()
    {
        X.clear();
        U.clear();
        Y.clear();
    }
    void get_iterator(BriefTrajectoryIterator<Model_> &traj_iter, size_t offset = 0)
    {
        traj_iter.Xiter = X.begin() + offset;
        traj_iter.Uiter = U.begin() + offset;
        traj_iter.Yiter = Y.begin() + offset;
    }
};

template <typename T, size_t xs, size_t us, size_t ys>
class Trajectory
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    Trajectory() {};
    Trajectory(T timeStep_, int horizon_)
    {        
        create_data(timeStep_, horizon_);
    }
    void create_data(T timeStep_, int horizon_);
    void zero_all();
    void zero_val_approx();
    void clear();
    void update_nominal_vals();

    void pop_front();
    void push_back_zero();
    int size(){return Xbar.size();}

public:
    T duration;
    T timeStep;
    int horizon; // the size of the state trajectory would be then horizon + 1

    /* Shared pointers to hold state, control and output trajectory */
    deque<VecM<T, xs>> Xbar;
    deque<VecM<T, xs>> X;
    deque<VecM<T, us>> Ubar;
    deque<VecM<T, us>> U;
    deque<VecM<T, ys>> Y;

    /* Shared pointers to hold linearized dynamics */
    deque<MatMN<T, xs, xs>> A;
    deque<MatMN<T, xs, us>> B;
    deque<MatMN<T, ys, xs>> C;
    deque<MatMN<T, ys, us>> D;

    /* Shared pointers to hold value function approximation */
    deque<T> V;
    deque<T> dV;
    deque<VecM<T, us>> dU;
    deque<VecM<T, xs>> G;
    deque<MatMN<T, xs, xs>> H;
    deque<MatMN<T, xs, xs>> K;

    deque<RCostData<T, xs, us, ys>> rcostData;
    TCostData<T, xs> tcostData;
};


template <typename EigenMat>
void set_eigen_deque_zero(deque<EigenMat> &mat_deque)
{
    for (auto &m : mat_deque)
    {
        m.setZero();
    }
}

template <typename T>
void set_scalar_deque_zero(deque<T> &scalar_deque)
{
    for (auto &v : scalar_deque)
    {
        v = 0;
    }
}

#endif