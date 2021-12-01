#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "utils/rk4.h"

namespace plant{

class Plant{
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
	double dlai_dt(PlantAssimilationResult& res, Env &env, PlantParameters &par, PlantTraits &traits);

	// demographics
	double p_survival_germination();
	
	double mortality_rate();

	double fecundity_rate(double mass, PlantTraits &traits);

	template<class Env>
	std::vector<double> growth_rates(Env &env);


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

			this->geometry->set_state(S.begin()+1, traits);
			prod = S[0];
			rep = S[4];
			seed_pool = S[5];
			germinated = S[6];

			std::vector<double> ss = growth_rates(env);
			dSdt[0] = ss[0]; //dmass_dt;	   // biomass production rate
			dSdt[1] = ss[1]; //dL_dt;       // lai growth rate
			dSdt[2] = ss[2]; //dsize_dt;    // size (diameter) growth rate
			dSdt[3] = ss[3]; //dlitter_dt;  // litter biomass growth rate
			dSdt[4] = ss[4]; //(1-fg)dBdt;  // litter biomass growth rate
			dSdt[5] = fecundity_rate(dSdt[4], traits) - seed_pool/par.ll_seed;
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
