#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "utils/rk4.h"

namespace plant{

class Plant{
	public:
	struct{
		// ** core demographic rates **
		double dlai_dt;
		double dsize_dt;
		double dmort_dt;
		double dseeds_dt;
		
		// intermediate and consistency-check variables
		double dmass_dt_tot;
		double dmass_dt_lai;
		double dmass_dt_rep;
		double dmass_dt_growth;
		double dmass_dt_lit;
	} rates;

	public:

	PlantTraits traits;
	PlantParameters par;

	Assimilator * assimilator;
	PlantGeometry * geometry;
	
	Plant();
	~Plant();

	int initParamsFromFile(std::string file);
	
	void set_size(double x);

	double get_biomass();
	
	// LAI model
	template<class Env>
	double dlai_dt(PlantAssimilationResult& res, Env &env);

	// demographics
	template<class Env>
	double p_survival_germination(Env &env);
	
	template<class Env>
	double mortality_rate(Env &env);

	template<class Env>
	double fecundity_rate(double mass, Env &env);

	template<class Env>
	void calc_growth_rates(Env &env);


	void print();

	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - grows plant over dt with constant assimilation rate A
	// ** 
	template<class Env>
	void grow_for_dt(double t, double dt, Env &env, double &prod, double &rep, double &seed_pool, double &germinated){

		auto derivs = [&env, &prod, &rep, &seed_pool, &germinated, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			//if (fabs(t - 2050) < 1e-5) 
			env.updateClimate(t);

			prod = S[0];
			this->geometry->set_state(S.begin()+1, traits);
			rep = S[4];
			seed_pool = S[5];
			germinated = S[6];

			calc_growth_rates(env);
			
			dSdt[0] = rates.dmass_dt_tot;	   // biomass production rate
			dSdt[1] = rates.dlai_dt;       // lai growth rate
			dSdt[2] = rates.dsize_dt;    // size (diameter) growth rate
			dSdt[3] = rates.dmass_dt_lit;  // litter biomass growth rate
			dSdt[4] = rates.dmass_dt_rep; //(1-fg)dBdt;  // reproduction biomass growth rate
			dSdt[5] = fecundity_rate(rates.dmass_dt_rep, traits) - seed_pool/par.ll_seed;
			dSdt[6] = seed_pool/par.ll_seed;
		};

		std::vector<double> S = {prod, geometry->lai, geometry->get_size(), geometry->litter_pool, rep, seed_pool, germinated};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		prod = S[0];
		geometry->set_state(S.begin()+1, traits);
		rep = S[4];
		seed_pool = S[5];
		germinated = S[6];

	}

			
};


}	// namespace plant

#include "../src/plant.tpp"

#endif
