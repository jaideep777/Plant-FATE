#include <iomanip>
#include <fstream>
#include <random>

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"

//#include "pspm_interface.h"

using namespace std;

class FixedEnvironment : public env::Climate, public env::LightEnvironment {
	public:
	void print(double t){
		Climate::print(t);
		LightEnvironment::print();
	}
};


class LifeHistoryOptimizer{
	public:
	plant::Plant P;
	FixedEnvironment C;

	double dt = 0.1; 

	double rep;
	double litter_pool;
	double seeds;
	double prod;

	public:

	void init(){
		rep = 0;
		litter_pool = 0;
		seeds = 0;
		prod = 0;

		C = FixedEnvironment();
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

	void printHeader(ostream &lfout){
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

	void printState(double t, ostream& lfout){
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

	void set_state(vector<double>::iterator it){
		P.geometry.set_lai(*it++);
		P.set_size(*it++);
		prod = *it++;
		litter_pool = *it++;
		rep = *it++;
//		P.state.seed_pool = *it++;
		seeds = *it++;
		P.state.mortality = *it++;
	}

	void get_rates(vector<double>::iterator it){
		*it++ = P.rates.dlai_dt;       // lai growth rate
		*it++ = P.rates.dsize_dt;    // size (diameter) growth rate
		*it++ = P.bp.dmass_dt_tot;	   // biomass production rate
		*it++ = P.bp.dmass_dt_lit;  // litter biomass growth rate
		*it++ = P.bp.dmass_dt_rep; //(1-fg)dBdt;  // reproduction biomass growth rate
//		*it++ = P.rates.dseeds_dt_pool;
		*it++ = P.rates.dseeds_dt;
		*it++ = P.rates.dmort_dt;
	}

	void grow_for_dt(double t, double dt){

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


	// void lho_set_traits(std::vector<double> tvec){
	// 	P.traits.c = tvec[0];
	// 	P.traits.a = tvec[1];
	// 	P.coordinateTraits();
	// }

	// std::vector<double> lho_get_traits(){
	// 	vector<double> tvec({
	// 		P.traits.c,
	// 		P.traits.a
	// 	});
	// 	return tvec;
	// }

	double calcFitness(){
		init();
		// lho_set_traits(tvec);
		for (double t=2000; t<=2500; t=t+dt){
			grow_for_dt(t, dt);
		}
		return seeds;
	}

};





int main(){


// Creating random number generator for soil water potential
	// ofstream fswp("swp.txt");

	// double prng_mean = -3.0;
	// double prng_stddev = -4.0;
	// std::default_random_engine generator;
	// std::normal_distribution<double> dist(prng_mean, prng_stddev);
	// for (double t=2000; t<=2100; t=t+10){
	// 	C.t_met.push_back(t);
	// 	double val = dist(generator);
	// 	if(val>0) val = 0;
	// 	C.v_met_swp.push_back(val);
	// 	fswp << C.v_met_swp[t] << "\n";
	// }
	
	LifeHistoryOptimizer lho;
	lho.init();
	lho.C.print(0);

	ofstream fout("assim1.txt");
	double dt = 0.1;
	lho.printHeader(fout);
	for (double t=2000; t<=2500; t=t+dt){
		lho.grow_for_dt(t, dt);
		lho.printState(t+dt, fout);
	}
	fout.close();
	// fswp.close();

	cout << "At last t: " << "\n" 
		 << setprecision(12) 
		 << "  Total biomass    = " << lho.P.get_biomass() << "\n"
		 << "  Total litter     = " << lho.litter_pool << "\n"
		 << "  Total reproduc   = " << lho.rep << "\n"
		 << "  Total bio+lit+rep = " << lho.P.get_biomass() + lho.litter_pool + lho.rep << "\n"
		 << "  Total production = " << lho.prod << "\n";
	

	double rel_error = abs((lho.P.get_biomass()+lho.litter_pool+lho.rep)/lho.prod - 1);
	cout << "Relative error in biomass accounting = " << rel_error << endl;
	if (rel_error > 2e-5) return 1;
	
	double fitness = lho.calcFitness();
	cout << "After calcFitness(): " << "\n" 
		 << setprecision(12) 
		 << "  Total biomass    = " << lho.P.get_biomass() << "\n"
		 << "  Total litter     = " << lho.litter_pool << "\n"
		 << "  Total reproduc   = " << lho.rep << "\n"
		 << "  Total bio+lit+rep = " << lho.P.get_biomass() + lho.litter_pool + lho.rep << "\n"
		 << "  Total production = " << lho.prod << "\n";
	cout << "Fitness = " << fitness << endl;


	return 0;
}

//# R script to test:
//dat = read.delim("/home/jaideep/codes/plant_fate_ppa/assim.txt")
//dat$heartwood_fraction = 1-dat$sapwood_fraction

//par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

//plot(dat$height~dat$diameter)
//plot(dat$crown_area~I(dat$height*dat$diameter))
//plot(dat$crown_area~dat$height)
//plot(dat$crown_area~dat$diameter)
//plot(dat$leaf_area~dat$crown_area)
//plot(dat$sapwood_fraction~dat$height)
//# plot(sqrt(4*dat$crown_area/pi)~dat$height)

//plot(dat$height~dat$i)
//plot(dat$diameter~dat$i)
//plot(dat$total_mass~dat$i)
//points(dat$total_prod~dat$i, type="l", col="red")

//matplot(y=cbind(dat$assim_gross,
//                dat$assim_net), 
//        x=dat$i, col=c("green3", "green4"), log="", pch=20)

//matplot(y=cbind(dat$rr,
//                dat$rs,
//                dat$rl), 
//        x=dat$i, col=c("yellow2", "yellow3", "yellow4"), log="", pch=20)

//matplot(y=cbind(dat$tr,
//                dat$tl), 
//        x=dat$i, col=c("orange3", "orange4"), log="", pch=20)


//D = dat$diameter
//H = dat$height
//-log(1-H/20)*20/D
