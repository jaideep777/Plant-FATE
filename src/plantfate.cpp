#include "plantfate.h"
#include <filesystem>
using namespace std;

Simulator::Simulator(std::string params_file) : I(params_file), S("IEBT", "rk45ck") {
	paramsFile = params_file; // = "tests/params/p.ini";
	I.readFile();

	parent_dir = I.get<string>("outDir");
	expt_dir   = I.get<string>("exptName");
	
	save_state = (I.get<string>("saveState") == "yes")? true : false;

	state_outfile  = I.get<string>("savedStateFile");
	config_outfile = I.get<string>("savedConfigFile");

	continueFrom_stateFile = I.get<string>("continueFromState");
	continueFrom_configFile = I.get<string>("continueFromConfig");
	continuePrevious = (continueFrom_configFile != "null") && (continueFrom_stateFile != "null");
	saveStateInterval = I.getScalar("saveStateInterval");

	traits_file = I.get<string>("traitsFile");
	n_species = I.getScalar("nSpecies");
	evolve_traits = (I.get<string>("evolveTraits") == "yes")? true : false;

	timestep = I.getScalar("timestep");  // ODE Solver timestep
 	delta_T = I.getScalar("delta_T");    // Cohort insertion timestep

	solver_method = I.get<string>("solver");
	res = I.getScalar("resolution");

	T_seed_rain_avg = I.getScalar("T_seed_rain_avg");

	T_return = I.getScalar("T_return");

	E.metFile = I.get<string>("metFile");
	E.co2File = I.get<string>("co2File");
	E.update_met = (E.metFile == "null")? false : true;
	E.update_co2 = (E.co2File == "null")? false : true;
	E.use_ppa = true;

	traits0.init(I);
	par0.init(I);
}


void Simulator::set_metFile(std::string metfile){
	E.metFile = metfile;
	E.update_met = (metfile == "")? false : true;
}


void Simulator::set_co2File(std::string co2file){
	E.co2File = co2file;
	E.update_co2 = (co2file == "")? false : true;
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

	y0 = tstart; //I.getScalar("year0");
	yf = tend;   //I.getScalar("yearf");
	ye = y0 + 120;  // year in which trait evolution starts (need to allow this period because r0 is averaged over previous time)

	// ~~~~~~~ Set up environment ~~~~~~~~~~~~~~~
	// E.metFile = met_file;
	// E.co2File = co2_file;
	E.init();
	E.print(0);

	// ~~~~~~~~~~ Create solver ~~~~~~~~~~~~~~~~~~~~~~~~~
	S = Solver(solver_method, "rk45ck");
	S.control.abm_n0 = 20;
	S.control.ode_ifmu_stepsize = timestep; //0.02; //0.0833333;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ifmu_order = 1;
	S.control.ebt_ucut = 1e-7;
	S.use_log_densities = true;
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
		//int nspp = I.getScalar("nSpecies");
		// int res = I.getScalar("resolution");
		for (int i=0; i<n_species; ++i){
			addSpeciesAndProbes(y0, 
								Tr.species[i].species_name, 
								Tr.species[i].lma, 
								Tr.species[i].wood_density, 
								Tr.species[i].hmat, 
								Tr.species[i].p50_xylem);
		}

		S.resetState(y0);
		S.initialize();
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


void Simulator::calc_seed_output(double t, Solver& S){
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
void Simulator::calc_r0(double t, double dt, Solver& S){
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		double r0 = log(spp->seeds_hist.get()/spp->birth_flux_in)/dt;
		
		spp->set_inputBirthFlux(spp->seeds_hist.get());
		spp->r0_hist.push(t, r0);
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

void Simulator::addSpeciesAndProbes(double t, string species_name, double lma, double wood_density, double hmat, double p50_xylem){

	PSPM_Plant p1;
	//p1.initFromFile(paramsFile);
	p1.par = par0;
	p1.traits = traits0;
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
	S.addSpecies(res, 0.01, 10, true, spp, 2, 1e-3);
	//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);

	// Add variants (probes) to solver
	if (evolve_traits){
		for (auto m : static_cast<MySpecies<PSPM_Plant>*>(spp)->probes) 
			S.addSpecies(res, 0.01, 10, true, m, 2, 1e-3);
	}

	S.copyCohortsToState();
}


void Simulator::simulate(){

	auto after_step = [this](double t){
		calc_seed_output(t, S);
		calc_r0(t, timestep, S);
	};

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
			if (t > ye){
				for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->calcFitnessGradient();
				for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->evolveTraits(delta_T);
			}
		}

		// Remove dead species
		vector<MySpecies<PSPM_Plant>*> toRemove;
		for (int k=0; k<S.species_vec.size(); ++k){
			auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
			if (spp->isResident){
				if (cwm.n_ind_vec[k] < 1e-6 && (t-spp->t_introduction) > 50) toRemove.push_back(spp);
			}
		}
		for (auto spp : toRemove) removeSpeciesAndProbes(spp);

		// // Shuffle species in the species vector -- just for debugging
		// if (int(t) % 10 == 0){
		// 	cout << "shuffling...\n";
		// 	std::random_shuffle(S.species_vec.begin(), S.species_vec.end());
		// 	S.copyCohortsToState();
		// }

		// Invasion by a random new species
		if (int(t) % 100 == 0){
			cout << "**** Invasion ****\n";
			addSpeciesAndProbes(t, 
			                    "spp_t"+to_string(t), 
			                    runif(0.05, 0.25),    //Tr.species[i].lma, 
			                    runif(300, 900),   //Tr.species[i].wood_density, 
			                    runif(2, 35),      //Tr.species[i].hmat, 
			                    runif(-6, -0.5)   //Tr.species[i].p50_xylem);
			);
		}

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
			double t_int = -log(double(rand())/RAND_MAX) * T_return;
			t_clear = t + fmin(t_int, 1000);
		}

		// Save simulation state at specified intervals
		if (int(t) % saveStateInterval == 0){
			saveState(&S, 
	          out_dir + "/" + std::to_string(t) + "_" + state_outfile, 
			  out_dir + "/" + std::to_string(t) + "_" + config_outfile, 
			  paramsFile);
		}

	}
	
}

