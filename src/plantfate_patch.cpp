#include "plantfate_patch.h"
#include <filesystem>
#include <csvrow.h>
using namespace std;

namespace pfate{

Patch::Patch(std::string params_file) : S("IEBT", "rk45ck") {
	config.paramsFile = params_file; // = "tests/params/p.ini";
	io::Initializer I;
	I.parse(params_file);

	config.parent_dir = I.get<string>("outDir");
	config.expt_dir   = I.get<string>("exptName");
	
	config.save_state = (I.get<string>("saveState") == "yes")? true : false;

	config.state_outfile  = I.get<string>("savedStateFile");
	config.config_outfile = I.get<string>("savedConfigFile");

	config.continueFrom_stateFile = I.get<string>("continueFromState");
	config.continueFrom_configFile = I.get<string>("continueFromConfig");
	config.continuePrevious = (config.continueFrom_configFile != "null") && (config.continueFrom_stateFile != "null");
	config.saveStateInterval = I.get<double>("saveStateInterval");

	config.traits_file = I.get<string>("traitsFile");
	config.n_species = I.get<double>("nSpecies");

	config.evolve_traits = (I.get<string>("evolveTraits") == "yes")? true : false;
	config.evolvable_traits = I.get_vector<std::string>("evolvableTraits");
	config.trait_scalars = I.get_vector<double>("traitScalars");
	config.trait_variances = I.get_vector<double>("traitVariances");
	config.T_r0_avg = I.get<double>("T_r0_averaging");
	if (config.trait_scalars.size() != config.evolvable_traits.size()) throw std::runtime_error("Please specify as many trait scalars as evolvable traits");
	if (config.trait_variances.size() != config.evolvable_traits.size()) throw std::runtime_error("Please specify as many trait variances as evolvable traits");

	config.timestep = I.get<double>("timestep");  // timestep
	config.time_unit = I.get_verbatim("time_unit");
	config.T_cohort_insertion = I.get<double>("T_cohort_insertion");    // Cohort insertion timestep

	config.solver_method = I.get<string>("solver");
	config.res = I.get<double>("resolution");

	config.T_invasion = I.get<double>("T_invasion");
	config.T_seed_rain_avg = I.get<double>("T_seed_rain_avg");
	config.T_return = I.get<double>("T_return");

	climate_stream.i_metFile = I.get<string>("i_metFile");
	climate_stream.a_metFile = I.get<string>("a_metFile");
	climate_stream.co2File   = I.get<string>("co2File");
	climate_stream.update_i_met = (climate_stream.i_metFile == "null")? false : true;
	climate_stream.update_a_met = (climate_stream.a_metFile == "null")? false : true;
	climate_stream.update_co2   = (climate_stream.co2File   == "null")? false : true;

	sio.emgProps_file = I.get<std::string>("emgProps");
	sio.cwmAvg_file = I.get<std::string>("cwmAvg");
	sio.cwmperSpecies_file = I.get<std::string>("cwmperSpecies");
	sio.traits_file = I.get<std::string>("traits");

	traits0.init(I);
	par0.init(I);
}


void Patch::set_i_metFile(std::string file){
	climate_stream.i_metFile = file;
	climate_stream.update_i_met = (file == "")? false : true;
}


void Patch::set_a_metFile(std::string file){
	climate_stream.a_metFile = file;
	climate_stream.update_a_met = (file == "")? false : true;
}


void Patch::set_co2File(std::string co2file){
	climate_stream.co2File = co2file;
	climate_stream.update_co2 = (co2file == "")? false : true;
}

std::vector<plant::PlantTraits> Patch::readTraitsFromFile(std::string fname){

	std::ifstream fin(fname.c_str());
	if (!fin){
		throw std::runtime_error("Could not open file " + fname + "\n");
	}
	
	std::vector<plant::PlantTraits> species;
	
	// read header
	flare::CSVRow row;
	fin >> row;
	
	plant::PlantTraits traits;
	// read time column for all rows
	while(fin >> row){
		std::string cell;

		cell = row[1];
		traits.species_name = cell;	

		cell = row[4];
		if (cell != "" && cell != "NA")	traits.wood_density = std::stod(cell)*1000; // g/cc to kg/m3
		else traits.wood_density = 686.638;

		cell = row[5];
		if (cell != "" && cell != "NA")	traits.hmat = std::stod(cell);
		else traits.hmat = 23.99;

		cell = row[6];
		if (cell != "" && cell != "NA")	traits.lma = std::stod(cell)*1e-3;  // convert g/m2 to kg/m2
		else traits.lma = 0.119378;

		cell = row[7];
		if (cell != "" && cell != "NA")	traits.p50_xylem = std::stod(cell); 
		else traits.p50_xylem = -2.29;

		// ignore further data (for now)	

		species.push_back(traits);		
	}

	return species;
}


void Patch::init(double tstart, double tend){
	config.out_dir  = config.parent_dir  + "/" + config.expt_dir; // TODO: Move to config::init()?

	// string command = "mkdir -p " + out_dir;
	std::filesystem::create_directories(config.out_dir);
	// string command2 = "cp " + paramsFile + " " + out_dir + "/p.ini";
	std::string copy_to = config.out_dir + "/p.ini";
	if (std::filesystem::exists(copy_to)) std::filesystem::remove(copy_to); // use this because the overwrite flag in below command does not work!
	std::filesystem::copy_file(config.paramsFile, copy_to, std::filesystem::copy_options::overwrite_existing);
	// int sysresult;
	// sysresult = system(command.c_str());
	// sysresult = system(command2.c_str());

	// ~~~~~~~ Set up time-points ~~~~~~~~~~~~~~~
	config.y0 = tstart; //I.get<double>("year0");
	config.yf = tend;   //I.get<double>("yearf");
	config.ye = config.y0 + config.T_r0_avg + 20;  // year in which trait evolution starts (need to allow this period because r0 is averaged over previous time)
	ts.set_units(config.time_unit);
	par0.set_tscale(ts.get_tscale());

	t_next_disturbance = config.T_return;
	t_next_invasion = config.T_invasion;
	t_next_savestate = config.y0; // this will write state once at the beginning too, which is probably unnecessary
	t_next_writestate = config.y0; // this will write state once at the beginning too, which is probably unnecessary

	// ~~~~~~~ Set up environment ~~~~~~~~~~~~~~~
	E.use_ppa = true;
	E.set_elevation(0);
	E.set_acclim_timescale(7); 
	climate_stream.init();

	// ~~~~~~~~~~ Create solver ~~~~~~~~~~~~~~~~~~~~~~~~~
	S = Solver(config.solver_method, "rk45ck");
	S.control.abm_n0 = 20;
	S.control.ode_ifmu_stepsize = 1e20; //timestep; //0.02; //0.0833333;
	S.control.cohort_insertion_dt = config.T_cohort_insertion;
	S.control.sync_cohort_insertion = false;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ebt_ucut = 1e-7;
	S.control.cm_use_log_densities = true;
	// S.control.cohort_insertion_tol = 1e-6;
	S.setEnvironment(&E);

	// Add species
	if (config.continuePrevious){
		restoreState(&S, config.continueFrom_stateFile, config.continueFrom_configFile);
		config.y0 = S.current_time; // replace y0
	}
	else {
		// ~~~~~~~~~~ Read initial trait values ~~~~~~~~~~~~~~~~~~~~~~~~~
		std::vector<plant::PlantTraits> file_species = readTraitsFromFile(config.traits_file);
		for (auto& s : file_species){
			std::cout << s.species_name << ":  " << s.lma << "\t" << s.wood_density << "\t" << s.hmat << "\t" << s.p50_xylem << "\n";
		}

		// ~~~ Create initial resident species pool from traits file ~~~~
		for (int i=0; i<config.n_species; ++i){
			plant::PlantTraits traits = traits0; 
			traits.species_name = file_species[i].species_name;
			traits.lma          = file_species[i].lma;
			traits.wood_density = file_species[i].wood_density;
			traits.hmat         = file_species[i].hmat;
			traits.p50_xylem    = file_species[i].p50_xylem; // runif(-3.5,-0.5);

			addSpeciesAndProbes(config.y0, traits);
		}

		// S.resetState(y0);
		S.initialize(config.y0);
	} 

//	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());

	S.print();	

	sio.S = &S;
	sio.openStreams(config.out_dir);
}


void Patch::close(){
	//S.print();
	sio.closeStreams();

	saveState(&S, 
	          config.out_dir + "/" + config.state_outfile, 
			  config.out_dir + "/" + config.config_outfile, 
			  config.paramsFile);

	// free memory associated
	for (auto s : S.species_vec) delete static_cast<AdaptiveSpecies<PSPM_Plant>*>(s); 

}


double Patch::runif(double rmin, double rmax){
	double r = double(rand())/RAND_MAX; 
	return rmin + (rmax-rmin)*r;
}


void Patch::removeSpeciesAndProbes(AdaptiveSpecies<PSPM_Plant>* spp){
	// delete species probes and remove their pointers from solver
	for (auto p : spp->probes){ // probes vector is not modified in the loop, so we can use it directly to iterate
		delete p;               // this will delete the object, but the pointer p and its copy in the solver remain
		S.removeSpecies(p);    // this will remove its pointer from the solver
	}

	// delete the resident itself and remove its pointer from solver
	delete spp;
	S.removeSpecies(spp);

	// update state vector
	S.copyCohortsToState();
}


void Patch::addSpeciesAndProbes(double t, const plant::PlantTraits& traits){

	PSPM_Plant p1;
	//p1.initFromFile(paramsFile);

	p1.init(par0, traits);

	((plant::Plant*)&p1)->print();
	
	//p1.geometry.set_lai(p1.par.lai0); // these are automatically set by init_state() in pspm_interface
	p1.set_size({0.01});
	
	AdaptiveSpecies<PSPM_Plant>* spp = new AdaptiveSpecies<PSPM_Plant>(p1);
	spp->species_name = traits.species_name;
	spp->evolvable_traits = config.evolvable_traits;
	spp->trait_scalars = config.trait_scalars; //{800, 25};
	spp->trait_variance = config.trait_variances; // vector<double>(2, 0.01);
	spp->r0_hist.set_interval(config.T_r0_avg);
	spp->t_introduction = t;

	spp->seeds_hist.set_interval(config.T_seed_rain_avg);

	if (config.evolve_traits) spp->createVariants(p1);

	// Add resident species to solver
	S.addSpecies({static_cast<int>(config.res)}, {0.01}, {10}, {true}, spp, 2, 1e-3);
	//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);

	// Add variants (probes) to solver
	if (config.evolve_traits){
		for (auto m : static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp)->probes) 
			S.addSpecies({static_cast<int>(config.res)}, {0.01}, {10}, {true}, m, 2, 1e-3);
	}

	S.copyCohortsToState();
}


void Patch::shuffleSpecies(){
	// Shuffle species in the species vector -- just for debugging
	cout << "shuffling...\n";
	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
	S.copyCohortsToState();
}


void Patch::removeDeadSpecies(double t){
	// Remove dead species
	vector<AdaptiveSpecies<PSPM_Plant>*> toRemove;
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k]);
		if (spp->isResident){
			if (cwm.n_ind_vec[k] < 1e-6 && (t-spp->t_introduction) > 50) toRemove.push_back(spp);
		}
	}
	for (auto spp : toRemove) removeSpeciesAndProbes(spp);
}


void Patch::addRandomSpecies(double t){
	cout << "**** Invasion ****\n";

	plant::PlantTraits traits = traits0;
	traits.species_name = "spp_t"+to_string(t);
	traits.lma          = runif(0.05, 0.25);    //Tr.species[i].lma, 
	traits.wood_density = runif(300, 900);   //Tr.species[i].wood_density, 
	traits.hmat         = runif(2, 35);      //Tr.species[i].hmat, 
	traits.p50_xylem    = runif(-6, -0.5);   //Tr.species[i].p50_xylem);

	addSpeciesAndProbes(t, traits);
}

void Patch::evolveTraits(double t, double dt_evolution){
	for (auto spp : S.species_vec) static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp)->calcFitnessGradient();
	for (auto spp : S.species_vec) static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp)->evolveTraits(dt_evolution);
}

void Patch::disturbPatch(double t){
	for (auto spp : S.species_vec){
		for (int i=0; i<spp->xsize(); ++i){
			auto& p = (static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp))->getCohort(i);
			p.geometry.lai = p.par.lai0;
			double u_new = spp->getU(i) * 0 * double(rand())/RAND_MAX;
			spp->setU(i, u_new);
		}
		spp->setX(spp->xsize()-1, spp->xb);
	}
	S.copyCohortsToState();
}


// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
/// @brief     Calculate growth rates of all species and update seed input
/// @param t   Current time 
/// @param dt  timestep over which species growth rate is to be calculated
/// @param S   Solver
/// @ingroup   trait_evolution
/// @details   Species growth rate is defined from the seed perspective, i.e., 
///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species) 
void Patch::calc_seedrain_r0(double t){
	
	// calculate output seed rain at time t, S(t)
	vector<double> seeds = S.newborns_out(t);
	for (int k=0; k<seeds.size(); ++k){
		if (seeds[k] < 0){
			cout << "seeds[" << k << "] = " << seeds[k] << endl;
			S.print();
			static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->seeds_hist.print_summary(); 
			cout.flush();
		}
	}

	// calculate r0 and update input seed rain
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k]);

		double dt = t - spp->seeds_hist.get_last_t(); // Note: Due to this line, first value of r0 will be garbage, unless initialized!
		if (dt < config.timestep/2){
			cout << "dt = " << dt << '\n';
			spp->r0_hist.print();
			spp->seeds_hist.print();
			throw std::runtime_error("r0 dt is less than the timestep!");
		}

		spp->seeds_hist.push(t, seeds[k]);        // Push S(t) into averager so that S_avg(t) can be computed

		double seeds_in = spp->birth_flux_in;     // This was S_avg(t-dt)
		double seeds_out = spp->seeds_hist.get(); // This is  S_avg(t),    i.e. average over seed-rain-avg interval, either a year or successional time window 
		double eps = 1e-20;
		double r0 = log((seeds_out+eps)/(seeds_in+eps))/dt;   // t0 = (log(S_avg(t)/S_avg(t-dt))/dt
		
		spp->set_inputBirthFlux(seeds_out);       // S_avg(t) will be input for next step
		spp->r0_hist.push(t, r0);                 // r0 is averaged again for better evolutoinary convergence 
		// spp->r0_hist.print_summary();
	}
}


/// @brief This function simulates patch to time t
/// @param t Final time up to which patch should be simulated
void Patch::simulate_to(double t){
	if (int(t*12) % 12 == 0){
		cout << "stepping = " << setprecision(6) << S.current_time << " --> " << t << "\t(";
		for (auto spp : S.species_vec) cout << spp->xsize() << ", ";
		cout << ")" << endl;
	}
	// Step size to be used for evolutionary dynamics, 
	// since trait evolution is done after completing step_to call
	double dt_evol = t - S.current_time; 

	// perform step
	auto after_step = [this](double _t){
		// Calc r0, update seed rain
		// calc_seedrain_r0(_t);
	};

	S.step_to(t, after_step);

	// update r0 and seed rain after step
	// this implies that seed tain and r0 are not updated during internal steps - I think thats okay.
	calc_seedrain_r0(t);

	// update output metrics - needed before removeDeadSpecies()
	cwm.update(t, S);
	props.update(t, S);

	// evolve traits
	if (config.evolve_traits && t > config.ye){
		evolveTraits(t, dt_evol);
	}

	// remove species whose total abundance has fallen below threshold (its probes are also removed)
	removeDeadSpecies(t); // needs updated cwm for species abundances

	// Invasion by a random new species
	if (t >= t_next_invasion){
		addRandomSpecies(t);
		t_next_invasion = t + config.T_invasion;
	}

	// clear patch by disturbance	
	if (t >= t_next_disturbance){
		disturbPatch(t);
		double dt_next = -log(double(rand())/RAND_MAX) * config.T_return;
		dt_next = std::clamp(dt_next, 0.0, 10*config.T_return);
		t_next_disturbance = t + dt_next;
	}

	// Save simulation state at specified intervals
	if (t > t_next_savestate || fabs(t-t_next_savestate) < 1e-6){
		saveState(&S, 
			config.out_dir + "/" + std::to_string(t) + "_" + config.state_outfile, 
			config.out_dir + "/" + std::to_string(t) + "_" + config.config_outfile, 
			config.paramsFile);
		
		t_next_savestate += config.saveStateInterval;
	}

	// Shuffle species - just for debugging. result shouldnt change
	// if (int(t) % 10 == 0) shuffleSpecies(); 
}


void Patch::update_climate(double julian_time, env::ClimateStream& c_stream){
	c_stream.updateClimate(julian_time, E);
}


void Patch::update_climate(double _co2, double _tc, double _vpd, double _ppfd, double _swp){
	E.clim_inst.tc   = _tc;
	E.clim_inst.ppfd = _ppfd;
	E.clim_inst.rn   = _ppfd/2;  // Tentative, placeholder
	E.clim_inst.vpd  = _vpd;
	E.clim_inst.co2  = _co2;
	E.clim_inst.swp  = _swp;
}

void Patch::update_climate_acclim(double t_julian, double _co2, double _tc, double _vpd, double _ppfd, double _swp){
	env::Clim C = E.clim_acclim;
	C.tc   = _tc;
	C.ppfd = _ppfd;
	C.rn   = _ppfd/2;  // Tentative, placeholder
	C.vpd  = _vpd;
	C.co2  = _co2;
	C.swp  = _swp;
	E.set_forcing_acclim(t_julian, C);
}


/// @brief Simulate patch dynamics
/// TODO: Evntually, this function should be moved to a higher controller, which can simulate multiple patches and manage data IO
/// In Plant-FATE, we should always simulate step-by-step because we are explicitly managing seed rain feedback
// To simulate step-by-step 
// --------------------------
//   1) set the following solver properties:
//	     --> S.control.ode_ifmu_stepsize = 1e20; // this will ensure that simulate_to() takes only 1 internal step
//	     --> S.control.cohort_insertion_dt = T_cohort_insertion; // this will insert cohorts automatically at the specified interval
//	     --> S.control.sync_cohort_insertion = false;  // this will prevent cohort insertion at the end of simulate_to()
//   2) simulate in a loop with increments of `timestep`, e.g.,
//      for (double t=y0; t <= yf+1e-6; t=t+timestep) {
//      	updateClimate(S.current_time);
//          simulate_to(t);
//      }
//   3) Climate should be updated once at the beginning of the step, and not in computeEnv(). Consider the following sequence of events
//      Ideally, seed input S1 = seed production in interval t0-t1 ==> depends on avg {u,E,C} over t0-t1. Now since dt is small (~day), 
//      u,E dont change as much, but C can change substantially from t0 to t1. Hence, B(u0,E0,C)~B(u1,E1,C) so it can be evaluated at end of step in AfterStep(). 
//      however, since C changes abruptly at the end of the timestep, we should use C from the step beginnning, which is what plants see throughout 
//      the step. If C-update is put in computeEnv, newborns_out() will update C and we will get B(u1,E1,C1) instead of B(u1,E1,C0). Hence, 
//      better to manage C ourselves and update it once at the beginning of the timestep.
//
//         S0              S1              S2          ( v = seed rain from current vegetation that turns into seedlings)
//                                        ||                                                     ( |        )
//          -----.-------> || -----.----> |||||        ( . = seed ) --->  ( | = seedling ) --->  ( | = tree )
//        -v----/------------v----/-----------v---       
//        {..}_/           {...}_/          {......}   
//         t0              t1              t2
//        
//         climate:          C1           
//         C0               ---------------- C2
//         _________________                --------
//   4) This step-by-step approach is also better suited to spatial simulations, where the entire grid of 
//      climate input is read once before the beginning of the step. Otherwise, if each patch tried to read 
//      its own input, this will try to update entire grid whenever any patch wants to upfate (inefficient) 
//      or worse, create problems when patches are parallelized
//
// To simulate a long interval (we shouldnt use this in Plant-FATE), 
// ----------------------------
//   1) set the following solver properties:
//	     --> S.control.ode_ifmu_stepsize = timstep; // this will ensure that simulate_to() internally steps by one `timestep` at a time
//	     --> S.control.cohort_insertion_dt = T_cohort_insertion; // this will insert cohorts automatically at the specified interval
//	     --> S.control.sync_cohort_insertion = false;  // this will prevent cohort insertion at the end of simulate_to()
//   2) Simulate in one go: simulate_to(T_final) OR 
//      break up the simulation into any desired number of intervals, where interval is at least several times the solver internal `timestep` 
//      (this is ideal for accuracy and efficiency, but not strictly necessary). E.g.,
//      for (double t=y0; t <= yf+1e-6; t=t+T_long) {
//      	simulate_to(t);
//      }
void Patch::simulate(){

	for (double t=config.y0; t <= config.yf+1e-6; t=t+config.timestep) {  // 1e-6 ensures that last timestep to reach yf is actually executed
		// read forcing inputs
		// std::cout << "update Env (explicit)... t = " << S.current_time << ":\n";
		update_climate(ts.to_julian(S.current_time)+1e-6, climate_stream); // The 1e-6 is to ensure that when t coincides exactly with time in climate file, we ensure that the value in climate file is read by asking for a slightly higher t
		// ((env::Climate&)E).print(t);

		// simulate patch
		simulate_to(t);

		// write outputs
		if (t > t_next_writestate || fabs(t-t_next_writestate) < 1e-6){
			sio.writeState(t, cwm, props);
			t_next_writestate += 1;
		}

	}
}

} // namespace pfate