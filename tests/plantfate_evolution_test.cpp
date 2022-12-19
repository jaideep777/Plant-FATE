#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
using namespace std;

#include <solver.h>
#include "pspm_interface.h"
#include "trait_reader.h"
#include "community_properties.h"
#include "trait_evolution.h"


inline double runif(double rmin=0, double rmax=1){
	double r = double(rand())/RAND_MAX; 
	return rmin + (rmax-rmin)*r;
}


/// @brief     Calculate seed output of all species
/// @param t   Current time 
/// @param S   Solver
/// @ingroup   trait_evolution
/// @details   Species seed output rate is defined as,  
///            \f[S = \int_{x_b}^{x_m}{f(s)u(s)ds}\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species) 
void calc_seed_output(double t, Solver& S){
	vector<double> seeds = S.newborns_out(t);
	// cout << "t = " << fixed << setprecision(10) << t << ", Species r0s:\n";
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		spp->seeds_hist.push(t, seeds[k]);
		if (seeds[k] < 0){
			cout << "seeds[" << k << "] = " << seeds[k] << endl;
			S.print();
			spp->seeds_hist.print_summary(); cout.flush();
		}
		// cout << "   " << k << ": " <<  S.species_vec[k]->birth_flux_in << " --> " << seeds[k] << "/" << seeds_hist[k].get() << ", r0 = " << setprecision(8) << r0 << "\n";
	}
}


// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
/// @brief     Calculate growth rates of all species and update seed input
/// @param t   Current time 
/// @param dt  timestep over which growth rate is to be calculated
/// @param S   Solver
/// @ingroup   trait_evolution
/// @details   Species growth rate is defined from the seed perspective, i.e., 
///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species) 
void calc_r0(double t, double dt, Solver& S){
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		double r0 = log(spp->seeds_hist.get()/spp->birth_flux_in)/dt;
		
		spp->set_inputBirthFlux(spp->seeds_hist.get());
		spp->r0_hist.push(t, r0);
		// spp->r0_hist.print_summary();
	}
}

void removeSpeciesAndProbes(Solver* S, MySpecies<PSPM_Plant>* spp){
	// delete species probes and remove their pointers from solver
	for (auto p : spp->probes){ // probes vector is not modified in the loop, so we can use it directly to iterate
		delete p;               // this will delete the object, but the pointer p and its copy in the solver remain
		S->removeSpecies(p);    // this will remove its pointer from the solver
	}

	// delete the resident itself and remove its pointer from solver
	delete spp;
	S->removeSpecies(spp);

	// update state vector
	S->copyCohortsToState();
}

void addSpeciesAndProbes(Solver *S, string params_file, io::Initializer &I, double t, string species_name, double lma, double wood_density, double hmat, double p50_xylem){
	int res = I.getScalar("resolution");
	bool evolve_traits = (I.get<string>("evolveTraits") == "yes")? true : false;
	double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");

	PSPM_Plant p1;
	p1.initParamsFromFile(params_file);
	p1.traits.species_name = species_name;
	p1.traits.lma = lma;
	p1.traits.wood_density = wood_density;
	p1.traits.hmat = hmat;
	p1.traits.p50_xylem = p50_xylem; // runif(-3.5,-0.5);
	
	p1.coordinateTraits();

	((plant::Plant*)&p1)->print();
	
	//p1.geometry.set_lai(p1.par.lai0); // these are automatically set by init_state() in pspm_interface
	p1.set_size(0.01);
	
	MySpecies<PSPM_Plant>* spp = new MySpecies<PSPM_Plant>(p1);
	spp->species_name = species_name;
	spp->trait_scalars = {0.2, 700};
	spp->fg_dx = 0.01;
	spp->trait_variance = vector<double>(2, 0.1);
	spp->r0_hist.set_interval(100);
	spp->t_introduction = t;

	spp->seeds_hist.set_interval(T_seed_rain_avg);

	if (evolve_traits) spp->createVariants(p1);

	// Add resident species to solver
	S->addSpecies(res, 0.01, 10, true, spp, 2, 1e-3);
	//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);

	// Add variants (probes) to solver
	if (evolve_traits){
		for (auto m : static_cast<MySpecies<PSPM_Plant>*>(spp)->probes) 
			S->addSpecies(res, 0.01, 10, true, m, 2, 1e-3);
	}

	S->copyCohortsToState();
}


int main(){

	// ~~~~~~~~~~ Read Paramaters ~~~~~~~~~~~~~~~~~~~~~~~~~
	string paramsFile = "tests/params/p.ini";
	io::Initializer I(paramsFile);
	I.readFile();
	string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
	string command = "mkdir -p " + out_dir;
	string command2 = "cp tests/params/p.ini " + out_dir;
	int sysresult;
	sysresult = system(command.c_str());
	sysresult = system(command2.c_str());

	// ~~~~~~~~~~ Set up environment ~~~~~~~~~~~~~~~~~~~~~~~~~
	PSPM_Dynamic_Environment E;
	E.metFile = I.get<string>("metFile");
	E.co2File = I.get<string>("co2File");
	E.init();
	E.print(0);
	E.use_ppa = true;
	E.update_met = true;
	E.update_co2 = true;

	// ~~~~~~~~~~ Create solver ~~~~~~~~~~~~~~~~~~~~~~~~~
	string solver_method = I.get<string>("solver");
	Solver S(solver_method, "rk45ck");
	S.control.abm_n0 = 20;
    S.control.ode_ifmu_stepsize = I.getScalar("timestep"); //0.02; //0.0833333;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ifmu_order = 1;
	S.control.ebt_ucut = 1e-7;
    S.use_log_densities = true;
	S.setEnvironment(&E);
	
	// ~~~~~~~~~~ Read initial trait values ~~~~~~~~~~~~~~~~~~~~~~~~~
	TraitsReader Tr;
	Tr.readFromFile(I.get<string>("traitsFile"));
	Tr.print();

	// ~~~ Create initial resident species pool from traits file ~~~~
	int nspp = I.getScalar("nSpecies");
	// int res = I.getScalar("resolution");
	bool evolve_traits = (I.get<string>("evolveTraits") == "yes")? true : false;
	for (int i=0; i<nspp; ++i){
		addSpeciesAndProbes(&S, paramsFile, I,
		                    I.getScalar("year0"), 
		                    Tr.species[i].species_name, 
		                    Tr.species[i].lma, 
		                    Tr.species[i].wood_density, 
		                    Tr.species[i].hmat, 
		                    Tr.species[i].p50_xylem);
	}

	S.resetState(I.getScalar("year0"));
	S.initialize();

//	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());

	S.print();

	SolverIO sio;
	sio.S = &S;
	sio.openStreams(out_dir, I);


	// ~~~~~~~~~~ Set up seed rain calculation ~~~~~~~~~~~~~~~~~~~~~~~~~
	// double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
	// vector<MovingAverager> seeds_hist(S.species_vec.size());
	// for (auto& M : seeds_hist) M.set_interval(T_seed_rain_avg);

	double timestep = I.getScalar("timestep");
	auto after_step = [&S, timestep](double t){
		calc_seed_output(t, S);
		calc_r0(t, timestep, S);
	};

//	ofstream fout("fmu_PlantFATE.txt");

	SpeciesProps cwm;
	EmergentProps props; 

	// ~~~~~~~~~~ Simulate ~~~~~~~~~~~~~~~~~~~~~~~~~
	double t_clear = 105000;
	// t is years since 2000-01-01
	double y0, yf;
	y0 = I.getScalar("year0");
	yf = I.getScalar("yearf");
	double delta_T = I.getScalar("delta_T");
	for (double t=y0; t <= yf; t=t+delta_T) {
		cout << "stepping = " << setprecision(6) << S.current_time << " --> " << t << "\t(";
		for (auto spp : S.species_vec) cout << spp->xsize() << ", ";
		cout << ")" << endl;

		S.step_to(t, after_step);

		// debug: r0 calc can be done here, it should give approx identical result compared to when r0_calc is dont in preCompute
		// S.step_to(t); //, after_step);
		// if (t > y0) after_step(t);
		// if (t > y0) calc_r0(t, delta_T, S);
		// //S.print(); cout.flush();

		cwm.update(t, S);
		props.update(t, S);
			
		sio.writeState(t, cwm, props);
	
		// evolve traits
		if (evolve_traits){
			if (t > y0 + 120){
				for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->calcFitnessGradient();
				for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->evolveTraits(delta_T);
			}
		}

		// // Remove dead species
		// vector<MySpecies<PSPM_Plant>*> toRemove;
		// for (int k=0; k<S.species_vec.size(); ++k){
		// 	auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		// 	if (spp->isResident){
		// 		if (cwm.n_ind_vec[k] < 1e-6 && (t-spp->t_introduction) > 50) toRemove.push_back(spp);
		// 	}
		// }
		// for (auto spp : toRemove) removeSpeciesAndProbes(&S, spp);

		// // Shuffle species in the species vector -- just for debugging
		// if (int(t) % 10 == 0){
		// 	cout << "shuffling...\n";
		// 	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
		// 	S.copyCohortsToState();
		// }

		// // Invasion by a random new species
		// if (int(t) % 300 == 0){
		// 	cout << "**** Invasion ****\n";
		// 	addSpeciesAndProbes(&S, paramsFile, I,
		// 	                    t, 
		// 	                    "spp_t"+to_string(t), 
		// 	                    runif(0.05, 0.25),    //Tr.species[i].lma, 
		// 	                    runif(300, 900),   //Tr.species[i].wood_density, 
		// 	                    runif(2, 35),      //Tr.species[i].hmat, 
		// 	                    runif(-6, -0.5)   //Tr.species[i].p50_xylem);
		// 	);
		// }

		// clear patch after 50 year	
		if (t >= t_clear){
			for (auto spp : S.species_vec){
				for (int i=0; i<spp->xsize(); ++i){
					auto& p = (static_cast<MySpecies<PSPM_Plant>*>(spp))->getCohort(i);
					p.geometry.lai = p.par.lai0;
					double u_new = spp->getU(i) * 0 * double(rand())/RAND_MAX;
					spp->setU(i, u_new);
				}
				spp->setX(spp->xsize()-1, 0);
			}
			S.copyCohortsToState();
			double t_int = -log(double(rand())/RAND_MAX) * I.getScalar("T_return");
			t_clear = t + fmin(t_int, 1000);
		}
		
	}
	
	//S.print();
	sio.closeStreams();

	// free memory associated
	for (auto s : S.species_vec) delete static_cast<MySpecies<PSPM_Plant>*>(s); 
}

