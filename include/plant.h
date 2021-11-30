#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"

namespace plant{

class Plant{
	public:

	PlantTraits traits;
	PlantParameters par;

	Assimilator * assimilator;
	PlantGeometry * geometry;
	
	Plant(){
		assimilator = new Assimilator();
		geometry = new PlantGeometry();
	}	
	~Plant(){
		delete assimilator;
		delete geometry;
	}


	int initParamsFromFile(std::string file){
		int i = par.initFromFile(file);
		geometry->initGeometry(0.01, par, traits);
	}

	void set_size(double x){
		geometry->set_size(x, traits);
	}

	double get_biomass(){
		return geometry->total_mass(traits);
	}

	
	// LAI model
	template<class Env>
	double dlai_dt(PlantAssimilationResult& res, Env &env, PlantParameters &par, PlantTraits &traits){
		double lai_curr = geometry->lai;
		geometry->set_lai(lai_curr + par.dl);
		auto res_plus = assimilator->biomass_growth_rate(env, geometry, par, traits);
		geometry->set_lai(lai_curr);
		
		double dnpp_dL = (res_plus.npp - res.npp)/geometry->crown_area/par.dl;
		double dE_dL = (res_plus.trans - res.trans)/geometry->crown_area/par.dl;

		double dL_dt = par.response_intensity*(dnpp_dL - 0.001*dE_dL); //geometry->dlai_dt(traits);
		
		return dL_dt;
	}


	// demographics
	double p_survival_germination(){
	
	}

	double mortality_rate(){
	
	}

	double fecundity_rate(double mass, PlantTraits &traits){
		return mass/(4*traits.seed_mass);
	}


	template<class Env>
	std::vector<double> growth_rates(Env &env){
		
		auto res = assimilator->biomass_growth_rate(env, geometry, par, traits);	
		
		double dL_dt = dlai_dt(res, env, par, traits);

		double max_alloc_lai = par.max_alloc_lai*std::max(res.npp, 0.0); // if npp is negative, there can be no lai increment. if npp is positive, max 10% can be allocated to lai increment
		double dmass_dt_lai = geometry->dmass_dt_lai(dL_dt, max_alloc_lai, traits);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?
			
		double dmass_dt = std::max(res.npp, 0.0);  // total mass increment (geometric and lai-driven), 0 if npp is negative

		double dmass_dt_geom = dmass_dt - std::max(dmass_dt_lai, 0.0);	 // biomass change due to allometric growth. if LAI is decreasing, no mass increase due to LAI
		double dlitter_dt = std::max(-dmass_dt_lai, 0.0);	// biomass from leaf loss goes into litter. if LAI is decreasing, leaves lost go into litter
		
		double dmass_growth_dmass = 1-geometry->dreproduction_dmass(par, traits);
		double dsize_dt = geometry->dsize_dmass(traits) * dmass_growth_dmass * dmass_dt_geom; // size growth rate. NOTE: LAI increment is prioritized (already subtracted from npp above)
	
		return {dmass_dt, dL_dt, dsize_dt, dlitter_dt, (1-dmass_growth_dmass)*dmass_dt_geom};
	}


	void print(){
		std::cout << "Plant:\n";
		std::cout << "  height = " << geometry->height << "\n";
		std::cout << "  diameter = " << geometry->diameter << "\n";
		std::cout << "  crown_area = " << geometry->crown_area << "\n";
		std::cout << "  lai = " << geometry->lai << "\n";
	}

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


#endif
