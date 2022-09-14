#include "TrajectoryManagement.h"
#include<algorithm>

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::create_data(T timeStep_, int horizon_){  
    timeStep = timeStep_;
    horizon = horizon_;
    duration = timeStep * horizon;

    Xbar.assign(horizon+1, VecM<T, xs>::Zero());
    X.assign(horizon+1, VecM<T, xs>::Zero());
    Ubar.assign(horizon, VecM<T, us>::Zero());
    U.assign(horizon, VecM<T, us>::Zero());
    Y.assign(horizon, VecM<T, ys>::Zero());

    A.assign(horizon+1, MatMN<T, xs, xs>::Zero());
    B.assign(horizon, MatMN<T, xs, us>::Zero());
    C.assign(horizon, MatMN<T, ys, xs>::Zero());
    D.assign(horizon, MatMN<T, ys, us>::Zero());

    V.assign(horizon+1, 0);
    dV.assign(horizon+1, 0);
    dU.assign(horizon, VecM<T, us>::Zero());
    G.assign(horizon+1, VecM<T, xs>::Zero());
    H.assign(horizon+1, MatMN<T, xs, xs>::Zero());
    K.assign(horizon+1, MatMN<T, xs, xs>::Zero());

    rcostData.assign(horizon, RCostData<T,xs,us,ys>());
}

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::clear(){
    Xbar.clear();
    X.clear();
    Ubar.clear();
    U.clear();
    Y.clear();
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    V.clear();
    dV.clear();
    dU.clear();
    G.clear();
    H.clear();
    K.clear();
    rcostData.clear();

    horizon = 0;
    timeStep = 0;
    duration = 0;
}

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::zero_all(){
    set_eigen_deque_zero(Xbar);
    set_eigen_deque_zero(X);
    set_eigen_deque_zero(Ubar);
    set_eigen_deque_zero(U);
    set_eigen_deque_zero(Y);

    set_eigen_deque_zero(A);
    set_eigen_deque_zero(B);
    set_eigen_deque_zero(C);
    set_eigen_deque_zero(D);

    set_scalar_deque_zero(V);
    set_scalar_deque_zero(dV);
    set_eigen_deque_zero(dU);
    set_eigen_deque_zero(G);
    set_eigen_deque_zero(H);
    set_eigen_deque_zero(K);
}

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::zero_val_approx(){
    set_scalar_deque_zero(V);
    set_scalar_deque_zero(dV);
    set_eigen_deque_zero(dU);
    set_eigen_deque_zero(G);
    set_eigen_deque_zero(H);
    set_eigen_deque_zero(K);
}

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::update_nominal_vals(){
    std::copy(X.begin(), X.end(), Xbar.begin());
    std::copy(U.begin(), U.end(), Ubar.begin());
    // set_eigen_deque_zero(K);
}


template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::pop_front(){    
    Xbar.pop_front();
    X.pop_front();
    Ubar.pop_front();
    U.pop_front();
    Y.pop_front();

    A.pop_front();
    B.pop_front();
    C.pop_front();
    D.pop_front();

    V.pop_front();
    dV.pop_front();
    dU.pop_front();
    G.pop_front();
    H.pop_front();
    K.pop_front();

    rcostData.pop_front();
    horizon--;
}

template <typename T, size_t xs, size_t us, size_t ys>
void Trajectory<T,xs,us,ys>::push_back_zero(){
    
    Xbar.push_back(VecM<T, xs>::Zero());
    X.push_back(VecM<T, xs>::Zero());
    Ubar.push_back(VecM<T, us>::Zero());
    U.push_back(VecM<T,us>::Zero());
    Y.push_back(VecM<T,ys>::Zero());

    A.push_back(MatMN<T, xs, xs>::Zero());
    B.push_back(MatMN<T, xs, us>::Zero());
    C.push_back(MatMN<T, ys, xs>::Zero());
    D.push_back(MatMN<T, ys, us>::Zero());

    V.push_back(T(0));
    dV.push_back(T(0));
    dU.push_back(VecM<T, us>::Zero());
    G.push_back(VecM<T, xs>::Zero());
    H.push_back(MatMN<T, xs, xs>::Zero());
    K.push_back(MatMN<T, xs, xs>::Zero());

    rcostData.push_back(RCostData<T,xs,us,ys>());
    horizon++;
}


// template class TrajectorySequence<ModelInfo<double,24,24,0>,ModelInfo<double,24,24,0>>;
template class Trajectory<double,24,24,0>;
