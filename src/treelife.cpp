#include "treelife.h"
using namespace std;

void ErgodicEnvironment::print(double t){
	Climate::print(t);
	LightEnvironment::print();
}


void LifeHistoryOptimizer::init(){
	rep = 0;
	litter_pool = 0;
	seeds = 0;
	prod = 0;

	C.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	C.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	//C.init();
	C.z_star = {15, 10, 5, 0};
	C.canopy_openness = {1, exp(-0.5*1.8), exp(-0.5*3.5), exp(-0.5*5.5)};
	C.n_layers = C.z_star.size()-1;

	// We are tracking the life-cycle of a seed: how many seeds does a single seed produce (having gone through dispersal, germination, and plant life stages)
	P = plant::Plant();
	P.initParamsFromFile("tests/params/p.ini");
	P.geometry.set_lai(1);
	P.set_size(0.01);
	// Simulation below starts at seedling stage. So account for survival until seedling stage
	P.state.mortality = -log(P.p_survival_dispersal(C)*P.p_survival_germination(C)); // p{fresh seed is still alive after germination} = p{it survives dispersal}*p{it survives germination}

	double total_prod = P.get_biomass();
	cout << "Starting biomass = " << total_prod << "\n";
	cout << "Mortality until seedling stage = " << P.state.mortality << "\n";

}

void LifeHistoryOptimizer::printHeader(ostream &lfout){
	lfout << "i" << "\t"
		<< "ppfd" <<"\t"
		<< "assim_net" << "\t"
		<< "assim_gross" << "\t"
		<< "rl" << "\t"
		<< "rr" << "\t"
		<< "rs" << "\t"
		<< "tl" << "\t"
		<< "tr" << "\t"
		<< "dpsi" << "\t"
		<< "vcmax" << "\t"
		<< "transpiration" << "\t"
		<< "height" << "\t"	
		<< "diameter" << "\t"	
		<< "crown_area" << "\t"	
		<< "lai" << "\t"	
		<< "sapwood_fraction" << "\t"
		<< "leaf_mass" << "\t"
		<< "root_mass" << "\t"
		<< "stem_mass" << "\t"
		<< "coarse_root_mass" << "\t"
		<< "total_mass" << "\t"
		<< "total_rep" << "\t"
		// << "seed_pool" << "\t"
		// << "germinated" << "\t"
		<< "fitness" << "\t"
		<< "total_prod" << "\t"
		<< "litter_mass" << "\t"
		<< "mortality" << "\t"
		<< "mortality_inst" << "\t"
		<< "leaf_lifespan" << "\t"
		<< "fineroot_lifespan" << "\t" 
		<< "\n";
}

void LifeHistoryOptimizer::printState(double t, ostream& lfout){
	lfout << t << "\t"
			<< C.clim.ppfd << "\t"
			<< P.assimilator.plant_assim.npp << "\t" 
			<< P.assimilator.plant_assim.gpp << "\t" 
			<< P.assimilator.plant_assim.rleaf << "\t"
			<< P.assimilator.plant_assim.rroot << "\t"
			<< P.assimilator.plant_assim.rstem << "\t"
			<< P.assimilator.plant_assim.tleaf << "\t"
			<< P.assimilator.plant_assim.troot << "\t"
			<< P.assimilator.plant_assim.dpsi_avg << "\t" 
			<< P.assimilator.plant_assim.vcmax_avg << "\t" 
			<< P.assimilator.plant_assim.trans << "\t" 
			<< P.geometry.height << "\t"	
			<< P.geometry.diameter << "\t"	
			<< P.geometry.crown_area << "\t"	
			<< P.geometry.lai << "\t"	
			<< P.geometry.sapwood_fraction << "\t"	
			<< P.geometry.leaf_mass(P.traits) << "\t"	
			<< P.geometry.root_mass(P.traits) << "\t"	
			<< P.geometry.stem_mass(P.traits) << "\t"	
			<< P.geometry.coarse_root_mass(P.traits) << "\t"	
			<< P.get_biomass() << "\t"
			<< rep << "\t"
		//  << P.state.seed_pool << "\t"
		//  << germinated << "\t"
			<< seeds << "\t"
			<< prod << "\t"
			<< litter_pool << "\t"
			<< P.state.mortality << "\t"
			<< P.rates.dmort_dt << "\t"
			<< 1/P.assimilator.kappa_l << "\t"
			<< 1/P.assimilator.kappa_r << "\t"
			"\n";
}

void LifeHistoryOptimizer::set_state(vector<double>::iterator it){
	P.geometry.set_lai(*it++);
	P.set_size(*it++);
	prod = *it++;
	litter_pool = *it++;
	rep = *it++;
//		P.state.seed_pool = *it++;
	seeds = *it++;
	P.state.mortality = *it++;
}

void LifeHistoryOptimizer::get_rates(vector<double>::iterator it){
	*it++ = P.rates.dlai_dt;       // lai growth rate
	*it++ = P.rates.dsize_dt;    // size (diameter) growth rate
	*it++ = P.bp.dmass_dt_tot;	   // biomass production rate
	*it++ = P.bp.dmass_dt_lit;  // litter biomass growth rate
	*it++ = P.bp.dmass_dt_rep; //(1-fg)dBdt;  // reproduction biomass growth rate
//		*it++ = P.rates.dseeds_dt_pool;
	*it++ = P.rates.dseeds_dt;
	*it++ = P.rates.dmort_dt;
}

void LifeHistoryOptimizer::grow_for_dt(double t, double dt){

	auto derivs = [this](double t, std::vector<double>&S, std::vector<double>&dSdt){
		//if (fabs(t - 2050) < 1e-5) 
		//env.updateClimate(t);
		set_state(S.begin());
		P.calc_demographic_rates(C, t);
		
		// Override Plant-FATE fecundity calculations 
		// We need to explicitly include plant mortality here for fitness calcs
		double fec = P.fecundity_rate(P.bp.dmass_dt_rep, C);
		P.rates.dseeds_dt =  fec * exp(-P.state.mortality);  // Fresh seeds produced = fecundity rate * p{plant is alive}
		// P.rates.dseeds_dt_germ =   P.state.seed_pool/P.par.ll_seed;   // seeds that leave seed pool proceed for germincation

		get_rates(dSdt.begin());
	};

	std::vector<double> S = {P.geometry.lai, P.geometry.get_size(), prod, litter_pool, rep, seeds, P.state.mortality};
	RK4(t, dt, S, derivs);
	//Euler(t, dt, S, derivs);
	set_state(S.begin());
}	


double LifeHistoryOptimizer::calcFitness(){
	init();
	// lho_set_traits(tvec);
	for (double t=2000; t<=2500; t=t+dt){
		grow_for_dt(t, dt);
	}
	return seeds;
}
