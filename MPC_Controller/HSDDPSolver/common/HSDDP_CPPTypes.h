#ifndef HSDDP_CPPTYPES_H
#define HSDDP_CPPTYPES_H

#define EIGEN_NO_DEBUG 1 // Prevent Eigen from asserting dimension mismatch for speedup
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#define CALC_DYN_AND_PAR 0
#define CALC_PARTIALS_ONLY 1
#define CALC_DYNAMICS_ONLY 2
#define N_TIMESTEPS_MAX 110


#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Cholesky> // Cholesky decomposition for positive definitness analysis

#define PI 3.141592653589793238

//  #define TIME_BENCHMARK

#define int_T long long int

template<typename T, size_t m, size_t n>
using MatMN = Eigen::Matrix<T, m, n>;

template<typename T, size_t m>
using VecM = Eigen::Matrix<T, m, 1>;

template<typename T>
using DMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// Rotation Matrix
template <typename T>
using RotMat = typename Eigen::Matrix<T, 3, 3>;

// 2x1 Vector
template <typename T>
using Vec2 = typename Eigen::Matrix<T, 2, 1>;

// 3x1 Vector
template <typename T>
using Vec3 = typename Eigen::Matrix<T, 3, 1>;

// 4x1 Vector
template <typename T>
using Vec4 = typename Eigen::Matrix<T, 4, 1>;

template <typename T>
using Mat3 = typename Eigen::Matrix<T, 3, 3>;

template <typename T>
using Mat4 = typename Eigen::Matrix<T, 4, 4>;

// 4x1 Vector
template <typename T>
using Quat = typename Eigen::Matrix<T, 4, 1>;

// Spatial Vector (6x1, all subspaces)
template <typename T>
using SVec = typename Eigen::Matrix<T, 6, 1>;

// Spatial Transform (6x6)
template <typename T>
using SXform = typename Eigen::Matrix<T, 6, 6>;


template<typename T>
using DVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using Chol =  Eigen::LDLT<DMat<T>>;

template <typename Scalar_, int xsize_, int usize_, int ysize_>
struct ModelInfo
{
    typedef Scalar_ Scalar;
    enum
    {
        state_dim = xsize_,
        contrl_dim = usize_,
        output_dim = ysize_
    };
};

#endif //HSDDP_CPPTYPES_H
