#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "utils/rk4.h"
#include "utils/moving_average.h"

namespace plant{

class Plant{
	public:
	// ** core state variables **
	struct{
		//double lai;         // these are in geometry 
		//double size;
		double mortality;     // cummulative mortality
		double seed_pool = 0; // seed pool size  // FIXME: This needs to be a species characteristic if nonlinearities come in, e.g. environmentally dependent germination
	} state;
	
	// ** core rates **
	struct{
		double dlai_dt;
		double dsize_dt;
		double dmort_dt;
		double dseeds_dt_pool;
		double dseeds_dt_germ;
		double rgr;
	} rates;	
		
	// results of biomass partitioning
	struct {
		double dmass_dt_lai;
		double dmass_dt_rep;
		double dmass_dt_growth;
		double dmass_dt_lit;
		double dmass_dt_tot;
	} bp;

	PlantAssimilationResult res;

	// seed output history
	MovingAverager seeds_hist;

	public:

	PlantTraits traits;
	PlantParameters par;

	Assimilator assimilator; // to use pointers here, need to apply rule of 5
	PlantGeometry geometry;
	
	public:
//	Plant();
//	~Plant();
//	Plant(const Plant &P);  // we need a copy constructor to correctly set geometry and assimilator pointers

	int initParamsFromFile(std::string file);
	int coordinateTraits();
	
	void set_size(double x);

	double get_biomass();
	
	// LAI model
	template<class Env>
	double lai_model(PlantAssimilationResult& res, double _dmass_dt_tot, Env &env);

	template<class Env>
	void partition_biomass(double dm_dt_tot, double dm_dt_lai, Env &env);

	// demographics
	template<class Env>
	double size_growth_rate(double _dmass_dt_growth, Env &env);

	template<class Env>
	double mortality_rate(Env &env);

	template<class Env>
	double fecundity_rate(double _dmass_dt_rep, Env &env);

	template<class Env>
	void calc_demographic_rates(Env &env);

	template<class Env>
	double p_survival_germination(Env &env);

	template<class Env>
	double p_survival_dispersal(Env &env);


	void print();

	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - grows plant over dt with constant assimilation rate A
	// ** 
	template<class Env>
	void grow_for_dt(double t, double dt, Env &env, double &prod, double &rep, double &litter_pool, double &germinated){

		auto derivs = [&env, &prod, &rep, &litter_pool, &germinated, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			//if (fabs(t - 2050) < 1e-5) 
			//env.updateClimate(t);

			geometry.set_lai(S[0]);
			set_size(S[1]);
//			this->geometry.set_state(S.begin(), traits);
			prod = S[2];
			litter_pool = S[3];
			rep = S[4];
			state.seed_pool = S[5];
			germinated = S[6];

			calc_demographic_rates(env);
			
			dSdt[0] = rates.dlai_dt;       // lai growth rate
			dSdt[1] = rates.dsize_dt;    // size (diameter) growth rate
			dSdt[2] = bp.dmass_dt_tot;	   // biomass production rate
			dSdt[3] = bp.dmass_dt_lit;  // litter biomass growth rate
			dSdt[4] = bp.dmass_dt_rep; //(1-fg)dBdt;  // reproduction biomass growth rate
			dSdt[5] = rates.dseeds_dt_pool;
			dSdt[6] = seeds_hist.get(); //rates.dseeds_dt_germ;
		};

		std::vector<double> S = {geometry.lai, geometry.get_size(), prod, litter_pool, rep, state.seed_pool, germinated};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		geometry.set_lai(S[0]);
		set_size(S[1]);
//		geometry.set_state(S.begin(), traits);
		prod = S[2];
		litter_pool = S[3];
		rep = S[4];
		state.seed_pool = S[5];
		germinated = S[6];
		seeds_hist.push(t+dt, rates.dseeds_dt_germ);
		//seeds_hist.print();

	}


};


}	// namespace plant

#include "../src/plant.tpp"

#endif
