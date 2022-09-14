#ifndef CASADIGEN_H
#define CASADIGEN_H

#define int_T long long int

#include "comp_foot_pos_casadi.h"
#include "comp_foot_jacob_1_casadi.h"
#include "comp_foot_jacob_2_casadi.h"
#include "comp_foot_jacob_3_casadi.h"
#include "comp_foot_jacob_4_casadi.h"
#include "hkinodyn_casadi.h"
#include "hkinodyn_par_casadi.h"
#include <vector>

/*
  @brief: Get the numerical evaluation of a CasadiGen function and the output sparsity pattern
  @params: 
          arg: T ptr to an array of pointers whose element points to an input variable
          res: T ptr to an array of pointers whose element points to an output variable
          max_sz_res: maximum size of output variables
*/
template<typename T>
void casadi_interface(std::vector<T *> ARG, std::vector<T *> RES, int max_sz_res,
                      int f(const T **, T **, int_T *, T *, int),
                      const int_T *f_sparse_out(int_T),
                      int f_work(int_T *, int_T *, int_T *, int_T *));

#endif //CASADIGEN_H