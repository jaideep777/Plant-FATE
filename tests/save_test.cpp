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
#include "state_restore.h"


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
	// spp->fg_dx = 0.01;
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
	string command2 = "cp " + paramsFile + " " + out_dir + "/p.ini";
	int sysresult;
	sysresult = system(command.c_str());
	sysresult = system(command2.c_str());

	bool save_state = (I.get<string>("saveState") == "yes")? true : false;
	string state_outfile  = out_dir + "/" + I.get<string>("savedStateFile");
	string config_outfile = out_dir + "/" + I.get<string>("savedConfigFile");

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

	// std::random_shuffle(S.species_vec.begin(), S.species_vec.end());


	SolverIO sio;
	sio.S = &S;
	sio.openStreams(out_dir, I);

	// saveState(&S, paramsFile);

	// ~~~~~~~~~~ Set up seed rain calculation ~~~~~~~~~~~~~~~~~~~~~~~~~
	// double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
	// vector<MovingAverager> seeds_hist(S.species_vec.size());
	// for (auto& M : seeds_hist) M.set_interval(T_seed_rain_avg);

	double timestep = I.getScalar("timestep");
	auto after_step = [&S, timestep](double t){
		calc_seed_output(t, S);
		calc_r0(t, timestep, S);
	};


	SpeciesProps cwm;
	EmergentProps props; 

	// ~~~~~~~~~~ Simulate ~~~~~~~~~~~~~~~~~~~~~~~~~
	double t_clear = 105000;
	// t is years since 2000-01-01
	double y0, yf;
	y0 = I.getScalar("year0");
	yf = I.getScalar("yearf");
	double delta_T = I.getScalar("delta_T");
	double t;
	for (t=y0; t <= 1200; t=t+delta_T) {
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

		// Shuffle species in the species vector -- just for debugging
		if (int(t) % 10 == 0){
			cout << "shuffling...\n";
			std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
			S.copyCohortsToState();
		}

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
	

	cout << " ---------------------- Solver before save ---------------------\n";
	S.print();


	saveState(&S, state_outfile, config_outfile, paramsFile);
	for (auto s : S.species_vec) delete static_cast<MySpecies<PSPM_Plant>*>(s); 
	S.species_vec.clear();
	S.print();

	S = Solver(solver_method, "rk45ck");
	S.control.abm_n0 = 20;
    S.control.ode_ifmu_stepsize = I.getScalar("timestep"); //0.02; //0.0833333;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ifmu_order = 1;
	S.control.ebt_ucut = 1e-7;
    S.use_log_densities = true;
	S.setEnvironment(&E);
//	sio.S = &S;

	restoreState(&S, state_outfile, config_outfile);
	cout << " ---------------------- Solver after restore ---------------------\n";
	S.print();


	for (; t <= yf; t=t+delta_T) {
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

		// Shuffle species in the species vector -- just for debugging
		if (int(t) % 10 == 0){
			cout << "shuffling...\n";
			std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
			S.copyCohortsToState();
		}

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


// PARAMS FILE FOR testing state save/restore
// 1. create a reference dataset by commenting out the save and restore lines
//    -> set exptName to something_ref
// 2. rerun by uncommenting those lines
//    -> set exptName to something
// Compare results in R

// # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// # Input parameters for Simulations
// # If you change the order of parameters below, you will get what you deserve
// # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// > STRINGS 
// # # > DIR
// # Directories for data output
// # homeDir			/home/chethana/codes/Flare/tests		# home dir - no spaces allowed
// traitsFile      tests/data/Amz_trait_filled_HD.csv
// metFile         tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv
// co2File         tests/data/CO2_ELE_AmzFACE2000_2100.csv

// outDir  		pspm_output		    # output dir name
// exptName 	 	ELE_HD_save_state_with_evol							# expt name

// emgProps        AmzFACE_D_PFATE_ELE_HD.txt
// cwmAvg          AmzFACE_Y_mean_PFATE_ELE_HD.txt       
// cwmperSpecies   AmzFACE_Y_PFATE_ELE_HD.txt
// traits          traits_ELE_HD.txt

// solver          IEBT

// evolveTraits    yes

// saveEndState    yes
// savedState      pf_saved_state.txt
// savedConfig     pf_saved_config.ini

// > SCALARS
// # ** 
// # ** Solver parameters
// # **
// resolution     5
// timestep       0.1
// delta_T        1

// # **
// # ** Simulation parameters
// # **
// nSpecies       3
// year0          1000
// yearf          1350
// nPatches       3    

// # **
// # ** Core traits (default values)  
// # **
// lma            0.122 # 0.1163 # leaf mass per leaf area [kg/m2]
// zeta           0.20			# fine root mass per leaf area [kg/m2] 0.7
// fcr            0.47         # coarse root mass per unit stem mass
// hmat           29.18		# maximum height [m]
// fhmat          0.8          # height at reproductive maturity as fraction of hmat
// seed_mass      3.8e-5	    # [kg]
// wood_density   690	        # [kg/m3]
// p50_xylem      -2.29        # Xylem P50 [MPa]  -4.515

// # p50_leaf     -1.5		# Leaf hydraulic capacity [MPa]
// K_leaf         0.5e-16		# Leaf conductance [m]  (Tuned to gs)
// K_xylem        4e-16		# Leaf conductance [m]
// b_leaf         1			# Shape parameter of leaf vulnerabilty curve [-] (Ref: Joshi et al 2022)
// b_xylem        1            # Shape parameter of xylem vulnerabilty curve [-] 

// # **
// # ** Phydro paramaters  
// # **
// kphio          0.087       # Quantum yield efficiency
// alpha          0.095       # Cost of maintaining photosynthetic capacity (Ref: Joshi et al 2022, removed outliers Helianthus and Glycine)
// gamma          1.052       # Cost of maintaining hydraulic pathway  (Ref: Joshi et al 2022, removed outliers Helianthus and Glycine)     


// # **
// # ** Allocation and geometric paramaters  
// # **
// m   1.5	 # crown shape paramaters
// n   2
// fg  0.1  # upper canopy gap fraction

// a   75     # height-diameter allometry 114
// c   6000   # crown area allometry
// b   0	   # bark allometry


// # ** LAI model
// optimize_lai          0  # 1
// Cc                    0.3  # Leaf construction costs per unit mass 
// Chyd                  0.00
// response_intensity    3  # speed of LAI response. This is calibrated to give ~3 months response lag
// lai_deriv_step     1e-4  # stepsize to calculate profit derivative wrt LAI
// max_alloc_lai		0.5	 # max fraction of npp that can be allocated to LAI increment
// lai0                  1.8  # initial LAI

// # **
// # ** Leaf Economics
// # **
// # For all values below, Ref: Wang et al 2021 Leaf economics)
// les_u             768   # [dimensionless]
// les_k1             24.5   # g biomass / mol CO2 (see cbio below)
// les_k2           0.0864   # (mol-CO2/day) / (umol-CO2/s)
// les_hT_dH         65.33e3  # J mol-1
// les_hT_c          26.35   # - 
// les_molar_R       8.31    # J mol-1 K-1
// les_cc            13    # dimensionless    (Ref: Colin)

// # ** 
// # ** Respiration and turnover 
// # **
// rd  0.015                 # ratio of leaf dark respiration rate to vcmax [-]  (Ref: Colin)
// # rl  
// rr  0.6275                  # 5.3   # fine-root respiration rate [yr-1]  (kg/kg/yr) kg biomass lost = rr * kg in roots
// rs  0.032                # 0.16  # sapwood respiration rate [yr-1]    (kg/kg/yr) kg biomass lost = rs * kg in sapwood

// ll   0.5   # leaf lifespan [yr]
// lr   1     # fine root lifespan [yr]

// cbio 2.45e-2    # kg biomass per mol CO2
// y    0.75		# yield factor accounting for growth respiration Ref: Educated guess, also used in other models)

// k_light	0.5		# Light extinction coefficient

// # ** 
// # ** Demographics
// # **
// a_f1 0.15    # Max fraction of biomass invested into reproduction
// a_f2 10      # steepness of switch to reproduction

// ll_seed  15   # seed lifespan in dormant pool (seed avg dormancy)


// # **
// # ** Dispersal and germination
// # **
// Sd            1e-5     # probability of survival during dispersal
// npp_Sghalf    0.5      # 1.5 required productivity for 50% probability of survival during germination

// # **
// # ** Mortality
// # ** 
// # mI      0.002    # 0.002 natural mortality rate
// # mD      0.008    # 0.001 mortality rate due to diameter
// # mD_e    0.2     # exponent in mortality rate due to diameter
// # mS      1      # 0.0  scaling of rgr in mortality rate due to carbon starvation
// # mS0     1e-4   # mort = -log(mS0 + mS*rgr)

// # mS0    1    # used for mort = mS0*exp(-mS*rgr)
// # mS     10 # 5000      

// # Parameters of mortality function
// # exp(c0 + clnD*log(D) + cD*D + cWD*(wd*wd-cWD0*cWD0) + par.cS0*exp(-par.cS*bp.dmass_dt_tot))
// # c0     -5
// # clnD   -0.3
// # cD      0.2
// # cWD    -1.48
// # cWD0    0.690
// # cS0     0   # 1
// # cS      10

// # r = c0 + cL*log(L) + clnD*log(D) + cD*(D) + cG*log(rgr*D) + cWD*(WD - cWD0)
// # M = 1/(1+exp(-r))
//  c0     -5   # a
//  cL      0 # 0.15   # b
//  clnD    0 # -1   # c
//  cD      0.005  # d
//  cG      -2 # -0.35   # e
//  cWD    -0.05
//  cWD0    650
//  cS0     0

//  cD0     0.06  # cD0*D^2
//  cD1     0.05 # cD1*exp(-D/0.01)

//  m_alpha  0.0598
//  m_beta   18.7159
//  m_gamma  0.0094

// # **
// # ** Disturbance
// # **
// T_return              100       # return interval of disturbance
// T_seed_rain_avg       0.099       # years over which seed rain is averaged


// > ARRAYS
// # > PARAM_SWEEPS
// # parameter sets to loop over. These override any of the parameters set above.
// # all of these MUST BE SET AFTER particle-system initialization. 
// # 	  else, particleSystem will initialize with 1st value in these vectors by default 
// # c_full = 0.1 0.12 0.14 0.16 0.19 0.22 0.26 0.3 0.35 0.41 0.48 0.56 0.65 0.77 0.89 1.05 1.22 1.43 1.67 1.96 2.29 2.68 3.13 3.66 4.28 5 5.85 6.84 8 9.36 10.95 12.8  -1
// # cS_full = 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 -1
// # c = 0.1 0.14 0.19 0.26 0.35 0.48 0.65 0.89 1.22 1.67 2.29 3.13 4.28 5.85 8 10.95 -1
// # c_offset = 0.12 0.16 0.22 0.3 0.41 0.56 0.77 1.05 1.43 1.96 2.68 3.66 5 6.84 9.36 12.8 -1
// # bvec		0.0002 0.002 0.02 0.2 -1


