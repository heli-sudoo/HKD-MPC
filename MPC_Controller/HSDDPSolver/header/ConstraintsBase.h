#ifndef CONSTRAINT_BASE_H
#define CONSTRAINT_BASE_H

#include "HSDDP_CPPTypes.h"
#include "HSDDP_CompoundTypes.h"
#include <deque>
#include <string>
#include <algorithm>
#include <memory>
#include <iostream>

using std::vector;
using std::deque;
using std::string;
using std::shared_ptr;

// Path constraint data structure (at one time step)
template<typename T, size_t xs, size_t us, size_t ys>
struct IneqConstrData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // probably not needed in this structure
	
    T g;
    VecM<T, xs> gx;
    VecM<T, us> gu;
    VecM<T, ys> gy;  
	MatMN<T, xs, xs> gxx;
	MatMN<T, us, us> guu;
	MatMN<T, ys, ys> gyy;

	IneqConstrData(){
        g = 0;
        gx.setZero();
        gu.setZero();
        gy.setZero();
		gxx.setZero();
		guu.setZero();
		gyy.setZero();
    }
	static IneqConstrData<T,xs,us,ys> Zero(){
		IneqConstrData<T,xs,us,ys> data;
		return data;
	}
};

// Terminal constraint data structure (at one time step)
template<typename T, size_t xs>
struct TConstrData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T h;  
    VecM<T, xs> hx;
	MatMN<T, xs, xs> hxx;

	TConstrData(){
		h = 0;
		hx.setZero();
		hxx.setZero();
	}    
};

// Augmented Lagrangian parameter structure
template<typename T>
struct AL_Param_Struct
{
    T lambda = 0;
    T sigma = 0;
	void update_penalty(T beta){
		sigma *= beta;
	}
	void update_Lagrange(T h){
		lambda += h*sigma;
	}
};

// Reduced barrier parameter structure
template<typename T>
struct  REB_Param_Struct
{
    T delta = 0.1;
    T delta_min = 0.01;
    T eps = 1;
	void update_relax(T beta){
		delta *= beta;
		delta = fmax(delta, delta_min);
	}
	void update_weight(T beta){
		eps *= beta;
	}
};

template<typename T, size_t xs_, size_t us_, size_t ys_>
class PathConstraintBase
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	typedef VecM<T, xs_> State;
    typedef VecM<T, us_> Contrl;
    typedef VecM<T, ys_> Output;

    enum {xs = xs_,
          ys = ys_,
          us = us_};

	int size;
	int len;
	string name;
	T max_violation = 0;
	deque<vector<IneqConstrData<T,xs_,us_,ys_>>> data;	
	deque<vector<REB_Param_Struct<T>>> params;
	REB_Param_Struct<T> param_init;

	PathConstraintBase(){
		size = 0;
		len = 0;		
	}
    PathConstraintBase(const string& name_){
        size = 0;
		len = 0;
        name = name_;
    }
	PathConstraintBase(int size_, int len_, const string& name_){
		size = size_;
		len = len_;
		name = name_;
	}    
	/*
		@brief:	Allocate memory for constraint data. Always need to do this after a PathConstraintBase object is created  
	*/
	void create_data(){		
		for (int i = 0; i < len; i++)
		{
			data.push_back(vector<IneqConstrData<T,xs_,us_,ys_>>(size));
		}		
	}		
	void clear_data(){
		data.clear();
	}
	void resize_data(){
		clear_data();
		create_data();
	}

	/*
		@brief:	Allocate memory for ReB params initialize them 
	*/
	void initialize_params(const REB_Param_Struct<T>& param_init_){  
		param_init = param_init_;      
		params.clear();
		for (int k = 0; k < len; k++)
		{
			params.push_back(vector<REB_Param_Struct<T>>(size, param_init_));
		}		
	}
	/*
		@brief:	Allocate memory for ReB params initialize them using default value
	*/
	void initialize_params(){
		REB_Param_Struct<T> param_init_;
		param_init_.delta = 0.01;
		param_init_.delta_min = 0.001;
		param_init_.eps = 1;
		initialize_params(param_init_);
	}
	void reset_params(){
		initialize_params(param_init);
	}
	void update_params(T thresh, T beta_relax, T beta_weight){
		for (int k = 0; k < len; k++)
		{//for each time step
			for (int i = 0; i < size; i++)
			{//for each constraint
				if (data[k][i].g > -thresh) // If constraint satisfied, do nothing
				{
					continue;
				}
				params[k][i].update_weight(beta_weight);
				params[k][i].update_relax(beta_relax);
			}
			
		}
		
	}
	void update_horizon_len(int len_){
        len = len_;
    }  

	void update_constraint_size(int size_){
		size = size_;
	}
	void update_max_violation(int k){
		if (k == 0){
			max_violation = 0;
		}			
		T max_violation_k = 0;
		for (auto& c: data[k])
		{
			max_violation_k = std::min(max_violation_k, c.g);
		}
		max_violation = std::min(max_violation, max_violation_k);
	}
	// call update_max_violation in compute_violation in the derived class
	virtual void compute_violation(const State&, const Contrl&, const Output&, int k) = 0;

public:
	void pop_front(){
		data.pop_front();
		params.pop_front();
		len --;
	}
	void push_back(){
		data.push_back(vector<IneqConstrData<T,xs_,us_,ys_>>(size));
		params.push_back(params.back());
		len ++;
	}
	void pop_front_n(int n){
		for (int i = 0; i < n; i++){
			pop_front();
		}		
	}
	void push_back_n(int n){
		for (int i = 0; i < n; i++){
			push_back();
		}
		
	}	
};

template<typename T, size_t xs_>
class TerminalConstraintBase
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef VecM<T, xs_> State;	
	
    enum {xs = xs_};
	int size;
	string name;
	vector<TConstrData<T, xs_>> data;	
	vector<AL_Param_Struct<T>> params;
	AL_Param_Struct<T> param_init;
	T max_violation;

	TerminalConstraintBase(){size=0;}
	TerminalConstraintBase(const string& name_){name = name_;}
	TerminalConstraintBase(int size_, const string& name_){
		size = size_;
		name = name_;
	}
	void create_data(){		
		data = vector<TConstrData<T,xs_>>(size);
	}
	void clear_data(){
		data.clear();
	}
	void resize_data(){
		data.resize(size);
	}
	void update_constraint_size(int size_){
		size = size_;
	}
	void initialize_params(){
		AL_Param_Struct<T> param_init_;
        param_init_.lambda = 0;
        param_init_.sigma = 5;
		initialize_params(param_init_);
	}
	void initialize_params(const AL_Param_Struct<T>& param_init_){     
		params.clear();   
		param_init = param_init_;
		params = vector<AL_Param_Struct<T>>(size, param_init_);
	}
	void reset_params(){
		// initialize_params(param_init);
		for (auto &param : params)
		{
			param.sigma = param_init.sigma;
			// param.lambda = 0;
		}
		
	}
	void update_params(T thresh, T beta){
		for (int i = 0; i < size; i++)
		{// for each constraint
			// if (fabs(data[i].h) < thresh) // if constraint satisfied, do nothing
			// {
			// 	continue;
			// }
			if (fabs(data[i].h) > 0.01) // if too large, increase penalty
			{
				params[i].update_penalty(beta);
			}else // if not too large, update Lagrange multiplier
			{
				params[i].update_Lagrange(data[i].h);
			}						
		}		
	}
	void update_max_violation(){
		max_violation = 0.0;
		for (auto &c:data)
		{
			max_violation = std::max(max_violation, fabs(c.h));
		}
		
	}
	// call uppdate_max_violation in compute_violation
	virtual void compute_violation(const State&) = 0;
};

template<typename T, size_t xs, size_t us, size_t ys>
class ConstraintContainer
{
public:
	typedef VecM<T, xs> State;
    typedef VecM<T, us> Contrl;
    typedef VecM<T, ys> Output;

	vector<shared_ptr<PathConstraintBase<T, xs, us, ys>>> pathConstraints;
	vector<shared_ptr<TerminalConstraintBase<T, xs>>> terminalConstraints;

	T max_pconstr;
	T max_tconstr;

public:

	void add_pathConstraint(shared_ptr<PathConstraintBase<T, xs, us, ys>>  pathConstraint_){		
		pathConstraints.push_back(pathConstraint_);
		pcontrs_size += pathConstraint_->size;
	}
	void add_terminalConstraint(shared_ptr<TerminalConstraintBase<T, xs>> terminalConstraint_){
		terminalConstraints.push_back(terminalConstraint_);
		tconstr_size += terminalConstraint_->size;
	}
	// Get the name and size information for all path constraints
	void get_pathConstraintsInfo(vector<string>& name, vector<int>& size_){
		for (auto constraint : pathConstraints)
		{
			name.push_back(constraint->name);
			size_.push_back(constraint->size);
		}
		
	}
	// Get the name and size information for all terminal constraints
	void get_terminalConstraintsInfo(vector<string>& name, vector<int>& size_){
		for (auto constraint : terminalConstraints){
			name.push_back(constraint->name);
			size_.push_back(constraint->size);
		}
	}	
	void compute_path_constraints(const State& x, const Contrl& u, const Output& y, int k){
		for (auto pathConstraint : pathConstraints){
			pathConstraint->compute_violation(x,u,y,k);
		}
	}
	void compute_terminal_constraints(const State& x){
		for (auto terminalConstraint : terminalConstraints){
			terminalConstraint->compute_violation(x);
		}
	}
	// May need better implementation for speed up
	void get_path_constraints(vector<IneqConstrData<T,xs,us,ys>>& pconstrs_data, int k){
		// clear the vector
		pconstrs_data.clear();
		for (auto pathConstraint : pathConstraints)
		{
			append_vector(pconstrs_data, pathConstraint->data[k]);
		}		
	}
	// May need better implementation for speed up	
	void get_terminal_constraints(vector<TConstrData<T,xs>>& tconstrs_data){
		// clear the vector
		tconstrs_data.clear();
		for (auto terminalConstraint:terminalConstraints)
		{
			append_vector(tconstrs_data, terminalConstraint->data);
		}		
	}	
	T get_max_pconstrs(){
		max_pconstr = 0;
		for (auto pathConstraint : pathConstraints)
		{
			max_pconstr = std::min(max_pconstr, pathConstraint->max_violation);
		}
		return max_pconstr;
	}
	T get_max_tconstrs(){
		max_tconstr = 0;
		for (auto terminalConstraint:terminalConstraints)
		{
			max_tconstr = std::max(max_tconstr, terminalConstraint->max_violation);
		}
		return max_tconstr;
	}
	// May need better implementation for speed up	
	void get_reb_params(vector<REB_Param_Struct<T>>& reb_params, int k){
		// empty the vector
		reb_params.clear();
		for (auto pathConstraint : pathConstraints)
		{
			append_vector(reb_params, pathConstraint->params[k]);
		}
	}
	// May need better implementation for speed up	
	void get_al_params(vector<AL_Param_Struct<T>>& al_params){
		// empty the vector
		al_params.clear();
		for (auto terminalConstraint:terminalConstraints)
		{
			append_vector(al_params, terminalConstraint->params);
		}
	}	
	void update_reb_params(HSDDP_OPTION option){
		// for each path contraint object, update its reb params over the entire horizon
		for (auto pathConstraint : pathConstraints)
		{
			pathConstraint->update_params(option.pconstr_thresh, option.update_relax, option.update_ReB);
		}
	}
	void update_al_params(HSDDP_OPTION option){			
		for (auto terminalConstraint:terminalConstraints)
		{
			terminalConstraint->update_params(option.tconstr_thresh, option.update_penalty);
		}	
	}
	void pop_front_n(int n){
		// for each path contraint object, remove n elements in the front from the constraint data and params
		for (auto pathConstraint : pathConstraints)
		{
			pathConstraint->pop_front_n(n);
		}
	}
	void push_back_n(int n){
		// for each path contraint object, add n elements in the end to the constraint data and params
		for (auto pathConstraint : pathConstraints)
		{
			pathConstraint->push_back_n(n);
		}
	}
	void reset_params(){
		for (auto terminalConstraint:terminalConstraints)
		{
			terminalConstraint->reset_params();
		}	
		for (auto pathConstraint : pathConstraints)
		{
			pathConstraint->reset_params();
		}
	}

public:
	int pcontrs_size = 0;
	int tconstr_size = 0;
};

#endif // CONSTRAINT_BASE_H