#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include <solver.h>
#include "pspm_environment.h"
#include "pspm_plant.h"

vector<double> fmu_create_grid(double xmin, double xmax, double dxmin = 1e-2, double dxmax=0.1, double multiplier=1.05){
	vector <double> xvec;
	double x = xmin, dx = dxmin;
	while(x<xmax){
		xvec.push_back(x);	
		dx = std::min(dx*multiplier, dxmax);
		x += dx;
	}
	xvec.push_back(xmax);
	return xvec; 
}

vector<double> my_log_seq(double x0, double xf, int N){
	vector<double> grid;
	for (int i=0; i<N; ++i) grid.push_back(exp(log(x0) + (double(i)/(N-1))*(log(xf)-log(x0))));
	return grid;
}

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

vector<double> generateDefaultCohortSchedule(double max_time){

	vector<double> tvec;

	const double multiplier=0.2, min_step_size=1e-5, max_step_size=2;
	
	assert(min_step_size > 0 && "The minimum step size must be greater than zero");
	
	double dt = 0.0, time = 0.0;
	tvec.push_back(time);
	while (time <= max_time) {
		dt = exp2(floor(log2(time * multiplier)));
		time += min(max(dt, min_step_size), max_step_size);
		tvec.push_back(time);
	}

	// Drop the last time; that's not going to be needed:
	if (tvec.size() >=1) 	// JAI: added to avoid overflow warning
		tvec.resize(tvec.size() - 1);

	return tvec;
}



class SolverIO{
	public:
	int nspecies;
	Solver * S;

	vector <vector<ofstream>> streams;

	void openStreams(vector<string> varnames){
		varnames.insert(varnames.begin(), "u");
		varnames.insert(varnames.begin(), "X");
		
		for (int s=0; s < S->species_vec.size(); ++s){
			auto spp = S->species_vec[s];
			vector<ofstream> spp_streams;
			
			for (int i=0; i<varnames.size(); ++i){
				stringstream sout;
				sout << "species_" << s << "_" << varnames[i] << ".txt";
				cout << sout.str() << endl;
				ofstream fout(sout.str().c_str());
				spp_streams.push_back(std::move(fout));
			}
			streams.push_back(std::move(spp_streams));
		}
	}

	void closeStreams(){
		for (int s=0; s<streams.size(); ++s){
			for (int j=0; j<streams[s].size(); ++j){
				streams[s][j].close();
			}
		}
	}

	void writeState(){
		for (int s=0; s < S->species_vec.size(); ++s){
			auto spp = (Species<PSPM_Plant>*)S->species_vec[s];

			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << S->current_time << "\t";

			for (int j=0; j<spp->xsize(); ++j){
				auto& C = spp->getCohort(j);
				streams[s][0] << C.x[0] << "\t";
				streams[s][1] << C.u << "\t";
				streams[s][2] << C.vars.mortality << "\t";
				streams[s][3] << C.viable_seeds << "\t";
				streams[s][4] << C.vars.area_heartwood << "\t";
				streams[s][5] << C.vars.mass_heartwood << "\t";

			}
			
			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << "\n";
		}
	}
};


int main(int argc, char ** argv){
	
	
	LightEnvironment env(1);	
	env.light_profile.print();	
	
	PSPM_Plant p1;

	PSPM_Plant p2;
	p2.lma = 0.2625;
	p2.initParameters();
	p2.vars.height = p2.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p2.vars.area_leaf = p2.par.area_leaf_0; 


	PSPM_Plant p3;
	p3.lma = 0.4625;
	p3.initParameters();
	p3.vars.height = p3.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p3.vars.area_leaf = p3.par.area_leaf_0; 

	
	cout << p1 << endl;
	//cout << p << endl;

    Species<PSPM_Plant> s1(p1);
    Species<PSPM_Plant> s2(p2);
    Species<PSPM_Plant> s3(p3);
    //M.p = M.seed = p;
	s1.print(); 
	
	//exit(1);

    Solver S(SOLVER_FMU, "rk45ck");
    //S.control.ode_ifmu_stepsize = 0.1;
	//S.control.ifmu_centered_grids = false; //true;
	S.control.cm_use_log_densities = true;
	S.control.ode_eps = 1e-4;
	S.control.integral_interpolate = false;
	S.setEnvironment(&env);
	//    S.createSizeStructuredVariables({"mort", "fec", "heart_area", "heart_mass"});
    
	double ip_seed_rain = 1;
	if (argc > 1) ip_seed_rain = stod(argv[1]);
	cout << "Simulating with input seed rain = " << ip_seed_rain << endl;
 
 	S.addSpecies({fmu_create_grid(p1.vars.height, 20)}, &s1, 4, ip_seed_rain);
	S.addSpecies({fmu_create_grid(p1.vars.height, 20)}, &s2, 4, ip_seed_rain);
	S.addSpecies({fmu_create_grid(p1.vars.height, 20)}, &s3, 4, ip_seed_rain);
	
#ifndef USE_INIT_DIST
	// ensure that only 1st cohort has finite density
	for (int k=0; k<S.species_vec.size(); ++k) 
		for (int i=1; i<S.species_vec[k]->xsize(); ++i) S.species_vec[k]->setU(i,0);
	S.copyCohortsToState();	
#endif

	S.initialize();
	S.print();
	
	double tmax = 105.32;
	if (argc > 2) tmax = stod(argv[2]);
	vector <double> times = generateDefaultCohortSchedule(tmax);
	for (auto t : times) cout << t << " "; cout << endl;

	
	SolverIO sio;
	sio.S = &S;
	sio.openStreams({"mort", "fec", "heart", "sap"});

	ofstream fli("light_profile_ind_plant.txt");
	ofstream fseed("seed_rains.txt");
	
	vector <vector<double>> seeds_out(S.species_vec.size());

	for (size_t i=0; i < times.size(); ++i){

		S.step_to(times[i]);		
		
		vector<double> seeds = S.newborns_out(times[i]);
		for (int s=0; s< S.species_vec.size(); ++s){
			double S_D = 0.25;
			seeds_out[s].push_back(seeds[s] * S_D * env.patch_age_density(times[i]));
		}

		cout << times[i] << " " << S.species_vec[0]->xsize() << " ";
		for (int i=0; i<S.n_species(); ++i) cout << seeds[i] << " ";
		cout << " | " << env.light_profile.npoints << " | " << Cohort<PSPM_Plant>::np << "\n";

		fseed << times[i] << "\t";
		for (int i=0; i<S.n_species(); ++i) fseed << seeds[i] << "\t";
		fseed << "\n";


		vector<double> xl = myseq(0, 20, 200);
		for (auto h : xl) fli << env.canopy_openness(h) << "\t";
		fli << endl;

		sio.writeState();

	}
	
	fli.close();
	fseed.close();
	sio.closeStreams();
	int nga=0, nma=0, nfa=0, npa=0;

//	for (auto spp : S.species_vec) {spp->fnEvals(); nga += spp->ng; nma += spp->nm; nfa += spp->nf; npa += spp->np;}
//	cout << "Number of calls to p/g/m/f functions: " << npa << " " << nga << " " << nma << " " << nfa << endl;
	cout << "Number of calls to p/g/m/f (static) : " << Cohort<PSPM_Plant>::np << " " << Cohort<PSPM_Plant>::ng << " " << Cohort<PSPM_Plant>::nm << " " << Cohort<PSPM_Plant>::nf << endl;
	cout << "Number of calls to derivs           : " << S.odeStepper.get_fn_evals() << endl;

	for (int s=0; s< S.n_species(); ++s){
		cout << "Seed rain for Species " << s << " (Lindh 18) = " << pn::integrate_trapezium(times, seeds_out[s]) << endl;
	}

//	for (int s=0; s< S.n_species(); ++s){
//		auto spp = S.species_vec[s];
//		vector <double> fec_vec;
//		fec_vec.reserve(spp->xsize());
//		for (int i=0; i<spp->xsize(); ++i){
//			auto C = ((Species<PSPM_Plant>*)spp)->getCohort(i);
//			double patch_age_density = env.patch_age_density(times[i]);
//			double S_D = 0.25;
//			double output_seeds = spp->birth_flux_in * S_D * patch_age_density * C.viable_seeds;
//			//cout << times[i] << " " << M.input_seed_rain << " " << S_D << " " << patch_age_density << " " << (*itf) << " | " << output_seeds << endl;
//			fec_vec.push_back(output_seeds);
//		}
//		cout << "Seed rain for Species " << s << " (Falster 17) = " << pn::integrate_trapezium(times, fec_vec) << endl;
//	}
	

}

