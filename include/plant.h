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

	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - grows plant over dt with constant assimilation rate A
	// ** 
	template<class Env>
	void grow_for_dt(double t, double dt, Env &env, double &prod){

		auto derivs = [&env, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			env.updateClimate(t);

			this->geometry->set_state(S.begin()+1, traits);
			
			double dmass_dt = this->assimilator->biomass_growth_rate(env, this->geometry, this->par, this->traits);
			
			double dL_dt = this->geometry->dlai_dt(traits);
			double dmass_dt_lai = this->geometry->dmass_dt_lai(dL_dt, traits);  // biomass change resulting from LAI change
			double dmass_dt_geom = dmass_dt - std::max(dmass_dt_lai, 0.0);	 // biomass change due to allometric growth
			double dlitter_dt = std::max(-dmass_dt_lai, 0.0);	// biomass from leaf loss goes into litter
			double dsize_dt = this->geometry->dsize_dmass(traits) * dmass_dt_geom; // size growth rate

			dSdt[0] = dmass_dt;	   // biomass production rate
			dSdt[1] = dL_dt;       // lai growth rate
			dSdt[2] = dsize_dt;    // size (diameter) growth rate
			dSdt[3] = dlitter_dt;  // litter biomass growth rate
		};

		std::vector<double> S = {prod, geometry->lai, geometry->get_size(), geometry->litter_pool};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		prod = S[0];
		geometry->set_state(S.begin()+1, traits);
	}

			
};


}	// namespace plant


#endif
