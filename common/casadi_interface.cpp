#include "casadi_interface.h"
#include <cstdio>
#include <assert.h>
template<typename T>
void casadi_interface(std::vector<T *> ARG, std::vector<T *> RES, int max_sz_res,
                      int f(const T **, T **, int_T *, T *, int),
                      const int_T *f_sparse_out(int_T),
                      int f_work(int_T *, int_T *, int_T *, int_T *))
{
    T **arg = nullptr;
    T **res = nullptr;
    arg = new T* [ARG.size()];
    res = new T* [RES.size()]; 
           
    int_T sz_arg, sz_res, n_arg, n_res, sz_iw, sz_w, *iw = nullptr;
    T *w = nullptr;

    // get the size of each input and output variables
    f_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

    assert(ARG.size()<=sz_arg);
    assert(RES.size()<=sz_res);

    for (size_t idx_arg = 0; idx_arg < ARG.size(); idx_arg++)
    {
        arg[idx_arg] = ARG[idx_arg];
    }

    for (size_t idx_res = 0; idx_res < RES.size(); idx_res++)
    {
        res[idx_res] = new T[max_sz_res];
    }    
       
    iw = new int_T[0];
    w = new T[0];

    // get function output
    f((const T **)arg, res, iw, w, 1);

    // get sparsity pattern
    const int_T *sppattern;
    const int_T *rowinfo, *colinfo;
    int nnz;

    int nrow, ncol;
    for (size_t idx_res = 0; idx_res < sz_res; idx_res++)
    {
        sppattern = f_sparse_out(idx_res); // get the sparsity pattern (pointer to const array) for the ith output
        nrow = sppattern[0];
        ncol = sppattern[1];                

        colinfo = sppattern + 2;
        rowinfo = colinfo + ncol + 1;
        nnz = colinfo[ncol];

        // copy data from res to RES
        int nzidx = 0;
        for (int colidx = 0; colidx < ncol; colidx++)
        {   
            while (nzidx < colinfo[colidx+1])
            {
                RES[idx_res][rowinfo[nzidx]+nrow*colidx] = res[idx_res][nzidx];
                nzidx++;
            }
            
        }
        
    }
    
    if(nullptr != iw) delete [] iw;
    if(nullptr !=w) delete [] w;

    for (int idx_res = 0; idx_res < RES.size(); idx_res++)
    {
        if(nullptr != res[idx_res]) delete [] res[idx_res];
    }
    delete [] arg;
    delete [] res;
}

template void casadi_interface<double>(std::vector<double *> ARG, std::vector<double *> RES, int max_sz_res,
                                       int f(const double **, double **, int_T *, double *, int),
                                       const int_T *f_sparse_out(int_T),
                                       int f_work(int_T *, int_T *, int_T *, int_T *));