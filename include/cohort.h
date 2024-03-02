#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cassert>
#include "io_utils.h"

template<class Ind>
class Cohort : public Ind {
	public:
	static int np, ng, nm, nf;   // number of evaluations of demographic functions

	// std::array<double,dim> x; // Moved to IndividualBase
	double u = -9.9e20; //std::numeric_limits<double>::quiet_NaN(); // Best to avoid nan's because they cant be read via filestreams
	int id = 0;	

	int group_id = 0;   // ID of the spatial group this cohort belongs to
	int group_size = 0; // size of the spatial group this cohort belongs to

	double birth_time = 0;
	bool remove = false;
	
	bool need_precompute = true;

	// Default constructor simply calls Individual's default ctor
	// NOTE: This means user's Ind class will need ctor with no arguments
	Cohort() : Ind() {
	}

	// Construct a cohort from Individual using copy constructor of Individual
	Cohort(const Ind& _ind) : Ind(_ind){
		// std::cout << "Constructing cohort form Individual:\n"; 
		// _ind.print();
		// std::cout << '\n';
		// Ind::print();
		// std::cout << '\n';
		// print();
		// std::cout << '\n';
	}

	void print_xu(std::ostream &out = std::cout){
		out << std::setw(6)  << std::setprecision(4) << birth_time;
		for (auto xx : Ind::x) out << std::setw(12) << std::setprecision(4) << xx; 
		out << std::setw(12) << std::setprecision(4) << u; 
	}

	void print(std::ostream &out = std::cout){
		print_xu(out);
		Ind::print(out);
	}

	void set_size(const std::vector<double>& _x){
		std::copy(_x.begin(), _x.end(), Ind::x.begin());
		need_precompute = true; // when size is updated, next rate calc will need precompute
		Ind::set_size(Ind::x);
	}

	
	//  These are defined here so that precompute trigger can be 
	//  checked before calling user-defined function 
	void preCompute(double t, void * _env){
		++np;
		// std::cout << "cohort precompute: "; print(); std::cout << "\n";
		Ind::preCompute(t,_env);	
		need_precompute = false;   // Once precompute is called, no need to further precompute until necessary
	}
	

	decltype(Ind::x) growthRate(double t, void * _env){
		++ng;
		if (need_precompute) preCompute(t,_env);
		// std::cout << "cohort growthRate(): "; print(); std::cout << "\n";
		return Ind::growthRate(t,_env);
	}
	

	double mortalityRate(double t, void * _env){
		++nm;
		if (need_precompute) preCompute(t,_env);
		// std::cout << "cohort mortRate(): "; print(); std::cout << "\n";
		return Ind::mortalityRate(t,_env);	
	}
	

	double birthRate(double t, void * _env){
		++nf;
		if (need_precompute) preCompute(t,_env);
		// std::cout << "cohort birthRate: "; print(); std::cout << "\n";
		return Ind::birthRate(t,_env);	
	}
	

	void save(std::ostream &fout, int n_extra_vars){
		// Save/restore individual first (for metadata)
		Ind::save(fout);

		// Then save cohort (for cohort state). This way, Individual metadata will be available when set_size() and set_state() are called in restore()
		fout << "Cohort<Ind>::v2" << "   ";
		fout << std::make_tuple(
		          id
		        , group_id
		        , group_size
		        , birth_time
		        , remove
		        , need_precompute); // we actually need not save need_precompute, because set_size() will always set it to 1 during restore
		fout << to_vector(Ind::x) << "   " << u << "   ";

		std::vector<double> cumm_vars(n_extra_vars);
		auto it = cumm_vars.begin();
		Ind::get_accumulators(it);
		fout << cumm_vars << '\n';
	}


	void restore(std::istream &fin, int n_extra_vars){
		Ind::restore(fin);

		std::string s; fin >> s; // discard version number
		assert(s == "Cohort<Ind>::v2");
		fin >> id
		    >> group_id
		    >> group_size
		    >> birth_time
		    >> remove
		    >> need_precompute;

		std::vector<double> _x;
		fin >> _x >> u;
		set_size(_x);

		std::vector<double> cumm_vars(n_extra_vars);
		fin >> cumm_vars;
		auto it = cumm_vars.begin();
		Ind::set_accumulators(it);
	}
	
	// FIXME other env dependent rates should also check for precompute
};

template<class Model>
int Cohort<Model>::np = 0;

template<class Model>
int Cohort<Model>::ng = 0;

template<class Model>
int Cohort<Model>::nm = 0;

template<class Model>
int Cohort<Model>::nf = 0;

#endif

