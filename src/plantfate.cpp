#include "plantfate.h"
#include <filesystem>
using namespace std;

Simulator::Simulator(std::string params_file) : S("IEBT", "rk45ck") {
	paramsFile = params_file; // = "tests/params/p.ini";
	I.parse(params_file);

	parent_dir = I.get<string>("outDir");
	expt_dir   = I.get<string>("exptName");
	
	save_state = (I.get<string>("saveState") == "yes")? true : false;

	state_outfile  = I.get<string>("savedStateFile");
	config_outfile = I.get<string>("savedConfigFile");

	continueFrom_stateFile = I.get<string>("continueFromState");
	continueFrom_configFile = I.get<string>("continueFromConfig");
	continuePrevious = (continueFrom_configFile != "null") && (continueFrom_stateFile != "null");
	saveStateInterval = I.get<double>("saveStateInterval");

	traits_file = I.get<string>("traitsFile");
	n_species = I.get<double>("nSpecies");
	evolve_traits = (I.get<string>("evolveTraits") == "yes")? true : false;

	timestep = I.get<double>("timestep");  // ODE Solver timestep
 	T_cohort_insertion = I.get<double>("T_cohort_insertion");    // Cohort insertion timestep

	solver_method = I.get<string>("solver");
	res = I.get<double>("resolution");

	T_invasion = I.get<double>("T_invasion");
	T_seed_rain_avg = I.get<double>("T_seed_rain_avg");
	T_return = I.get<double>("T_return");

	climate_stream.metFile = I.get<string>("metFile");
	climate_stream.co2File = I.get<string>("co2File");
	climate_stream.update_met = (climate_stream.metFile == "null")? false : true;
	climate_stream.update_co2 = (climate_stream.co2File == "null")? false : true;

	E.use_ppa = true;

	traits0.init(I);
	par0.init(I);
}


void Simulator::set_metFile(std::string metfile){
	climate_stream.metFile = metfile;
	climate_stream.update_met = (metfile == "")? false : true;
}


void Simulator::set_co2File(std::string co2file){
	climate_stream.co2File = co2file;
	climate_stream.update_co2 = (co2file == "")? false : true;
}


void Simulator::init(double tstart, double tend){
	out_dir  = parent_dir  + "/" + expt_dir;

	// string command = "mkdir -p " + out_dir;
	std::filesystem::create_directories(out_dir);
	// string command2 = "cp " + paramsFile + " " + out_dir + "/p.ini";
	std::string copy_to = out_dir + "/p.ini";
	if (std::filesystem::exists(copy_to)) std::filesystem::remove(copy_to); // use this because the overwrite flag in below command does not work!
	std::filesystem::copy_file(paramsFile, copy_to, std::filesystem::copy_options::overwrite_existing);
	// int sysresult;
	// sysresult = system(command.c_str());
	// sysresult = system(command2.c_str());

	y0 = tstart; //I.get<double>("year0");
	yf = tend;   //I.get<double>("yearf");
	ye = y0 + 120;  // year in which trait evolution starts (need to allow this period because r0 is averaged over previous time)

	t_next_disturbance = T_return;
	t_next_invasion = T_invasion;
	t_last_evolution = 1e20;
	t_next_savestate = y0; // this will write state once at the beginning too, which is probably unnecessary
	t_next_writestate = y0; // this will write state once at the beginning too, which is probably unnecessary

	// ~~~~~~~ Set up environment ~~~~~~~~~~~~~~~
	// E.metFile = met_file;
	// E.co2File = co2_file;
	climate_stream.init();

	// ~~~~~~~~~~ Create solver ~~~~~~~~~~~~~~~~~~~~~~~~~
	S = Solver(solver_method, "rk45ck");
	S.control.abm_n0 = 20;
	S.control.ode_ifmu_stepsize = 1e20; //timestep; //0.02; //0.0833333;
	S.control.cohort_insertion_dt = T_cohort_insertion;
	S.control.sync_cohort_insertion = false;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ebt_ucut = 1e-7;
	S.control.cm_use_log_densities = true;
	S.setEnvironment(&E);

	// Add species
	if (continuePrevious){
		restoreState(&S, continueFrom_stateFile, continueFrom_configFile);
		y0 = S.current_time; // replace y0
	}
	else {
		// ~~~~~~~~~~ Read initial trait values ~~~~~~~~~~~~~~~~~~~~~~~~~
		TraitsReader Tr;
		Tr.readFromFile(traits_file);
		Tr.print();

		// ~~~ Create initial resident species pool from traits file ~~~~
		//int nspp = I.get<double>("nSpecies");
		// int res = I.get<double>("resolution");
		for (int i=0; i<n_species; ++i){
			plant::PlantTraits traits = traits0; 
			traits.species_name = Tr.species[i].species_name;
			traits.lma          = Tr.species[i].lma;
			traits.wood_density = Tr.species[i].wood_density;
			traits.hmat         = Tr.species[i].hmat;
			traits.p50_xylem    = Tr.species[i].p50_xylem; // runif(-3.5,-0.5);

			addSpeciesAndProbes(y0, traits);
		}

		// S.resetState(y0);
		S.initialize(y0);
	} 

//	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());

	S.print();	

	sio.S = &S;
	sio.openStreams(out_dir, I);
}


void Simulator::close(){
	//S.print();
	sio.closeStreams();

	saveState(&S, 
	          out_dir + "/" + state_outfile, 
			  out_dir + "/" + config_outfile, 
			  paramsFile);

	// free memory associated
	for (auto s : S.species_vec) delete static_cast<MySpecies<PSPM_Plant>*>(s); 

}


double Simulator::runif(double rmin, double rmax){
	double r = double(rand())/RAND_MAX; 
	return rmin + (rmax-rmin)*r;
}


// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
/// @brief     Calculate growth rates of all species and update seed input
/// @param t   Current time 
/// @param dt  timestep over which species growth rate is to be calculated
/// @param S   Solver
/// @ingroup   trait_evolution
/// @details   Species growth rate is defined from the seed perspective, i.e., 
///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species) 
void Simulator::calc_seedrain_r0(double t){
	
	// calculate output seed rain at time t, S(t)
	vector<double> seeds = S.newborns_out(t);
	for (int k=0; k<seeds.size(); ++k){
		if (seeds[k] < 0){
			cout << "seeds[" << k << "] = " << seeds[k] << endl;
			S.print();
			static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k])->seeds_hist.print_summary(); 
			cout.flush();
		}
	}

	// calculate r0 and update input seed rain
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);

		double dt = t - spp->seeds_hist.get_last_t(); // Note: Due to this line, first value of r0 will be garbage, unless initialized!

		spp->seeds_hist.push(t, seeds[k]);        // Push S(t) into averager so that S_avg(t) can be computed

		double seeds_in = spp->birth_flux_in;     // This was S_avg(t-dt)
		double seeds_out = spp->seeds_hist.get(); // This is  S_avg(t),    i.e. average over seed-rain-avg interval, either a year or successional time window 
		double r0 = log(seeds_out/seeds_in)/dt;   // t0 = (log(S_avg(t)/S_avg(t-dt))/dt
		
		spp->set_inputBirthFlux(seeds_out);       // S_avg(t) will be input for next step
		spp->r0_hist.push(t, r0);                 // r0 is averaged again for better evolutoinary convergence 
		// spp->r0_hist.print_summary();
	}
}


void Simulator::removeSpeciesAndProbes(MySpecies<PSPM_Plant>* spp){
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


void Simulator::addSpeciesAndProbes(double t, const plant::PlantTraits& traits){

	PSPM_Plant p1;
	//p1.initFromFile(paramsFile);

	p1.init(par0, traits);

	((plant::Plant*)&p1)->print();
	
	//p1.geometry.set_lai(p1.par.lai0); // these are automatically set by init_state() in pspm_interface
	p1.set_size({0.01});
	
	MySpecies<PSPM_Plant>* spp = new MySpecies<PSPM_Plant>(p1);
	spp->species_name = traits.species_name;
	spp->trait_scalars = {0.2, 700};
	// spp->fg_dx = 0.01;
	spp->trait_variance = vector<double>(2, 0.1);
	spp->r0_hist.set_interval(100);
	spp->t_introduction = t;

	spp->seeds_hist.set_interval(T_seed_rain_avg);

	if (evolve_traits) spp->createVariants(p1);

	// Add resident species to solver
	S.addSpecies({static_cast<int>(res)}, {0.01}, {10}, {true}, spp, 2, 1e-3);
	//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);

	// Add variants (probes) to solver
	if (evolve_traits){
		for (auto m : static_cast<MySpecies<PSPM_Plant>*>(spp)->probes) 
			S.addSpecies({static_cast<int>(res)}, {0.01}, {10}, {true}, m, 2, 1e-3);
	}

	S.copyCohortsToState();
}


void Simulator::shuffleSpecies(){
	// Shuffle species in the species vector -- just for debugging
	cout << "shuffling...\n";
	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
	S.copyCohortsToState();
}


void Simulator::removeDeadSpecies(double t){
	// Remove dead species
	vector<MySpecies<PSPM_Plant>*> toRemove;
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		if (spp->isResident){
			if (cwm.n_ind_vec[k] < 1e-6 && (t-spp->t_introduction) > 50) toRemove.push_back(spp);
		}
	}
	for (auto spp : toRemove) removeSpeciesAndProbes(spp);
}


void Simulator::addRandomSpecies(double t){
	cout << "**** Invasion ****\n";

	plant::PlantTraits traits = traits0;
	traits.species_name = "spp_t"+to_string(t);
	traits.lma          = runif(0.05, 0.25);    //Tr.species[i].lma, 
	traits.wood_density = runif(300, 900);   //Tr.species[i].wood_density, 
	traits.hmat         = runif(2, 35);      //Tr.species[i].hmat, 
	traits.p50_xylem    = runif(-6, -0.5);   //Tr.species[i].p50_xylem);

	addSpeciesAndProbes(t, traits);
}

void Simulator::evolveTraits(double t, double dt_evolution){
	for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->calcFitnessGradient();
	for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->evolveTraits(dt_evolution);
}

void Simulator::disturbPatch(double t){
	for (auto spp : S.species_vec){
		for (int i=0; i<spp->xsize(); ++i){
			auto& p = (static_cast<MySpecies<PSPM_Plant>*>(spp))->getCohort(i);
			p.geometry.lai = p.par.lai0;
			double u_new = spp->getU(i) * 0 * double(rand())/RAND_MAX;
			spp->setU(i, u_new);
		}
		spp->setX(spp->xsize()-1, spp->xb);
	}
	S.copyCohortsToState();
}



/// @brief This function simulates patch to time t
/// @param t Final time up to which patch should be simulated
void Simulator::simulate_to(double t){
	cout << "stepping = " << setprecision(6) << S.current_time << " --> " << t << "\t(";
	for (auto spp : S.species_vec) cout << spp->xsize() << ", ";
	cout << ")" << endl;

	// Step size to be used for evolutionary dynamics, 
	// since trait evolution is done after completing step_to call
	double dt_evol = t - S.current_time; 

	// perform step
	auto after_step = [this](double _t){
		// Calc r0, update seed rain
		calc_seedrain_r0(_t);
	};

	S.step_to(t, after_step);

	// update output metrics 
	cwm.update(t, S);
	props.update(t, S);

	// evolve traits
	if (evolve_traits && t > ye){
		evolveTraits(t, dt_evol);
	}

	// remove species whose total abundance has fallen below threshold (its probes are also removed)
	removeDeadSpecies(t); // needs updated cwm for species abundances

	// Invasion by a random new species
	if (t >= t_next_invasion){
		addRandomSpecies(t);
		t_next_invasion = t + T_invasion;
	}

	// clear patch by disturbance	
	if (t >= t_next_disturbance){
		disturbPatch(t);
		double dt_next = -log(double(rand())/RAND_MAX) * T_return;
		dt_next = std::clamp(dt_next, 0.0, 10*T_return);
		t_next_disturbance = t + dt_next;
	}

	// Save simulation state at specified intervals
	if (t >= t_next_savestate){
		saveState(&S, 
			out_dir + "/" + std::to_string(t) + "_" + state_outfile, 
			out_dir + "/" + std::to_string(t) + "_" + config_outfile, 
			paramsFile);
		
		t_next_savestate += saveStateInterval;
	}

	// Shuffle species - just for debugging. result shouldnt change
	// if (int(t) % 10 == 0) shuffleSpecies(); 
}


void Simulator::update_climate(double t, ClimateStream& c_stream){
	c_stream.updateClimate(flare::yearsCE_to_julian(t), E.clim);
}


void Simulator::update_climate(double t, double _co2, double _tc, double _vpd, double _ppfd, double _ppfd_max, double _swp){
	E.clim.tc = _tc;
	E.clim.ppfd_max = _ppfd_max;
	E.clim.ppfd = _ppfd;
	E.clim.vpd = _vpd;
	E.clim.co2 = _co2;
	E.clim.elv = _elv;
	E.clim.swp = _swp;
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
void Simulator::simulate(){

	for (double t=y0; t <= yf+1e-6; t=t+timestep) {
		// read forcing inputs
		climate_stream.updateClimate(flare::yearsCE_to_julian(S.current_time), E.clim);
		std::cout << "update Env (explicit)... t = " << S.current_time << ": tc = " << E.clim.tc << '\n';

		// simulate patch
		simulate_to(t);

		// write outputs
		if (t >= t_next_writestate){
			sio.writeState(t, cwm, props);
			t_next_writestate += 1;
		}

	}
}
