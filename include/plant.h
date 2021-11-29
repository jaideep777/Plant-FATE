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
		
			auto res = this->assimilator->biomass_growth_rate(env, this->geometry, this->par, this->traits);	
			
			double dmass_dt, dL_dt, dmass_dt_lai;
		   
			if (res.npp < 0){
				dmass_dt = dL_dt = dmass_dt_lai = 0;
			}
			else{
				dmass_dt = res.npp;

				double dl = 1e-4;
				double lai_back = this->geometry->lai;
				this->geometry->set_lai(lai_back+dl);
				auto res_plus = this->assimilator->biomass_growth_rate(env, this->geometry, this->par, this->traits);
				this->geometry->set_lai(lai_back);
				
				double dnpp_dL = (res_plus.npp - res.npp)/this->geometry->crown_area/dl;
				double dE_dL = (res_plus.trans - res.trans)/this->geometry->crown_area/dl;

				double response_intensity = 5;
				
				dL_dt = response_intensity*(dnpp_dL - 0.001*dE_dL); //this->geometry->dlai_dt(traits);
				dmass_dt_lai = std::min(this->geometry->dmass_dt_lai(dL_dt, traits), 0.1*dmass_dt);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?
				
				dL_dt = this->geometry->dlai_dt(dmass_dt_lai, traits);
			}
			
			double dmass_dt_geom = dmass_dt - dmass_dt_lai;	 // biomass change due to allometric growth
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
