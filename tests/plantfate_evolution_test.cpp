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


void calc_r0(double t, Solver& S, vector<MovingAverager>& seeds_hist){
	for (int k=0; k<S.species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
		double r0 = seeds_hist[k].get()/S.species_vec[k]->birth_flux_in;
		
		spp->set_inputBirthFlux(seeds_hist[k].get());
		spp->r0_hist.push(t, r0);
		// spp->r0_hist.print_summary();
	}
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
	E.update_co2 = false;

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

	// ~~~~~~~~~~ Create resident species ~~~~~~~~~~~~~~~~~~~~~~~~~
	int nspp = I.getScalar("nSpecies");
	int res = I.getScalar("resolution");
	// bool evolveTraits = (I.get<string>("evolveTraits") == "yes")? true : false;
	for (int i=0; i<nspp; ++i){
		PSPM_Plant p1;
		p1.initParamsFromFile("tests/params/p.ini");
		p1.traits.species_name = Tr.species[i].species_name;
		p1.traits.lma = Tr.species[i].lma;
		p1.traits.wood_density = Tr.species[i].wood_density;
		p1.traits.hmat = Tr.species[i].hmat;
		p1.traits.p50_xylem = Tr.species[i].p50_xylem; // runif(-3.5,-0.5);
		
		p1.coordinateTraits();

		((plant::Plant*)&p1)->print();
		
		//p1.geometry.set_lai(p1.par.lai0); // these are automatically set by init_state() in pspm_interface
		p1.set_size(0.01);
		
		MySpecies<PSPM_Plant>* spp = new MySpecies<PSPM_Plant>(p1);
		spp->trait_scalars = {0.2, 700};
		spp->fg_dx = 0.01;
		spp->trait_variance = vector<double>(2, 0.1);
		spp->r0_hist.set_interval(100);

		spp->createVariants(p1);

		S.addSpecies(res, 0.01, 10, true, spp, 3, 1e-3);
		//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);
		
		//	S.addSpecies(vector<double>(1, p1.geometry.get_size()), &spp, 3, 1);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	}

	// ~~~~~~~~~~ Create variant probes ~~~~~~~~~~~~~~~~~~~~~~~~~
	vector<Species_Base*> species_vec_copy = S.species_vec;
	for (auto spp : species_vec_copy){
		for (auto m : static_cast<MySpecies<PSPM_Plant>*>(spp)->probes) 
			S.addSpecies(res, 0.01, 10, true, m, 3, 1e-3);
	}

	S.resetState(I.getScalar("year0"));
	S.initialize();

	S.print();

	SolverIO sio;
	sio.S = &S;
	sio.openStreams(out_dir, I);


	// ~~~~~~~~~~ Set up seed rain calculation ~~~~~~~~~~~~~~~~~~~~~~~~~
	double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
	vector<MovingAverager> seeds_hist(S.species_vec.size());
	for (auto& M : seeds_hist) M.set_interval(T_seed_rain_avg);

	auto after_step = [&S, &seeds_hist](double t){
		vector<double> seeds = S.newborns_out(t);
		// cout << "t = " << fixed << setprecision(10) << t << ", Species r0s:\n";
		for (int k=0; k<S.species_vec.size(); ++k){
			seeds_hist[k].push(t, seeds[k]);
			if (seeds[k] < 0){
				cout << "seeds[" << k << "] = " << seeds[k] << endl;
				S.print();
				seeds_hist[k].print_summary(); cout.flush();
			}
			// cout << "   " << k << ": " <<  S.species_vec[k]->birth_flux_in << " --> " << seeds[k] << "/" << seeds_hist[k].get() << ", r0 = " << setprecision(8) << r0 << "\n";
		}

		calc_r0(t, S, seeds_hist);
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
		// if (t > y0) after_step(t);
		// if (t > y0) calc_r0(t, S, seeds_hist);
		//S.print(); cout.flush();

		cwm.update(t, S);
		props.update(t, S);
			
		sio.writeState(t, seeds_hist, cwm, props);
	
		if (t > y0 + 120){
			for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->calcFitnessGradient();
			for (auto spp : S.species_vec) static_cast<MySpecies<PSPM_Plant>*>(spp)->evolveTraits(delta_T);
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
			double t_int = -log(double(rand())/RAND_MAX) * I.getScalar("T_return");;
			t_clear = t + fmin(t_int, 1000);
		}
		
	}
	
	//S.print();
	sio.closeStreams();

	// free memory associated
	for (auto s : S.species_vec) delete static_cast<MySpecies<PSPM_Plant>*>(s); 
}

