#ifndef  PSPM_PSPM_INDIVIDUAL_BASE_H_
#define  PSPM_PSPM_INDIVIDUAL_BASE_H_

#include <iostream>
#include <vector>
#include <array>
#include <string>

#include "io_utils.h"

template<size_t dim = 1>
class IndividualBase{
	public:

	// Since we cant force users to manage x, these variables are managed by Cohort.
	// But it is defined here becuase it is templated and Cohort doesnt have access to `dim`
	std::array<double, dim> x; 
	// ---

	public:
	size_t state_size = dim;
	std::vector<std::string> varnames;
	
	// essential functions which must be defined by user
	virtual void set_size(const std::array <double, dim>& _x) = 0; // Note that Ind::set_size() is not guaranteed to set x.
	virtual double init_density(void * _env, double bf) = 0; 
	virtual std::array<double,dim> growthRate(double t, void * _env) = 0;
	virtual double mortalityRate(double t, void * _env) = 0;
	virtual double birthRate(double t, void * _env) = 0;

	// optional functions which can be defined by user through overloads
	virtual ~IndividualBase(){
	}
	
	virtual void preCompute(double t, void * _env){
	}

	virtual double establishmentProbability(double t, void  * _env){
		return 1;
	};

	// Functions related to additional state variables
	virtual void init_accumulators(double t, void * _env){}

	virtual std::vector<double>::iterator set_accumulators(std::vector<double>::iterator &it){
		return it;
	}
	virtual std::vector<double>::iterator get_accumulators(std::vector<double>::iterator &it){
		return it;
	}
	virtual std::vector<double>::iterator get_accumulatorRates(std::vector<double>::iterator &it){
		return it;
	} 

	// printing
	virtual void print(std::ostream &out = std::cout) const {}

	virtual void save(std::ostream& fout){}
	virtual void restore(std::istream& fin){}

};


#endif


