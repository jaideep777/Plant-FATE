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
		return par.initFromFile(file);
	}

	void set_size(double x){
		geometry->set_size(x, par, traits);
	}

	double get_biomass(){
		return geometry->total_mass(par, traits);
	}

	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - grows plant over dt with constant assimilation rate A
	// ** 
	template<class Env>
	void grow_for_dt(double t, double dt, Env &env, double &prod){

		auto derivs = [&env, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			this->geometry->set_size(S[1], par, traits);
			
			double dmass_dt = this->assimilator->biomass_growth_rate(1.0, env, this->geometry, this->par, this->traits);

			dSdt[0] = dmass_dt;	// biomass production rate
			dSdt[1] = this->geometry->dsize_dmass(par, traits) * dmass_dt; 
		};

		std::vector<double> S = {prod, geometry->get_size()};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		geometry->set_size(S[1], par, traits);
		prod = S[0];
	}

			
};


}	// namespace plant


#endif
