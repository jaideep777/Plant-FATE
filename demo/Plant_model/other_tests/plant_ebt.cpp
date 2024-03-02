#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
#include "solver.h"

#include "plant_pspm.h"

vector<double> generateDefaultCohortSchedule(double max_time){

	vector<double> tvec;

	const double multiplier=0.2, min_step_size=1e-5, max_step_size=2.0;
	
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


template<class Model, class Environment>
class SolverIO{
	public:
	int nspecies;
	Solver<Model, Environment> * S;

	vector <vector<ofstream>> streams;

	void openStreams(){
		for (int s=0; s < S->n_species(); ++s){
			auto spp = S->get_species(s);
			vector<string> varnames = spp->get_varnames();
			vector<ofstream> spp_streams;
			for (string name : varnames){
				stringstream sout;
				sout << "species_" << s << "_" << name + ".txt";
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
		for (int s=0; s < S->n_species(); ++s){
			auto spp = S->get_species(s);
			auto iset = spp->get_iterators(S->state);
			auto& ivec = iset.get();

			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << S->current_time << "\t";

			for (iset.rbegin(); !iset.rend(); --iset){
				for (int i=0; i<ivec.size(); ++i){
					streams[s][i] << *ivec[i] << "\t"; 
				}
			}
			
			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << "\n";
		}
	}
};


int main(){
	
	//initPlantParameters(plant::par);
	
	LightEnvironment env(1);	
	env.light_profile.print();	
	
	//plant::Environment env(1);
	plant::Plant p1;

	plant::Plant p;
	p.lma = 0.2625;
	p.initParameters();
	p.vars.height = p.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p.vars.area_leaf = p.par.area_leaf_0; 


	plant::Plant p3;
	p3.lma = 0.4625;
	p3.initParameters();
	p3.vars.height = p3.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p3.vars.area_leaf = p3.par.area_leaf_0; 

	
	cout << p1 << endl;
	cout << p << endl;

	//exit(1);

    Solver<PlantModel, LightEnvironment> S(SOLVER_EBT);
    S.control.cm_use_log_densities = true;
	S.control.ode_eps = 1e-4;
	S.setEnvironment(&env);
	//    S.createSizeStructuredVariables({"mort", "fec", "heart_area", "heart_mass"});
    
	PlantModel M1;
	M1.p = M1.seed = p1;
    cout << "HT1 === " << M1.p.vars.height << endl;
	

    PlantModel M;
    M.p = M.seed = p;
    cout << "HT === " << M.p.vars.height << endl;


    PlantModel M3;
    M3.p = M3.seed = p3;
    cout << "HT === " << M3.p.vars.height << endl;

	S.addSpecies(vector<double>(2, M1.p.vars.height), &M1, {"mort", "fec", "heart", "sap"}, M1.input_seed_rain);
	S.addSpecies(vector<double>(2, M.p.vars.height), &M, {"mort", "fec", "heart", "sap"}, M.input_seed_rain);
	S.addSpecies(vector<double>(2, M3.p.vars.height), &M3, {"mort", "fec", "heart", "sap"}, M3.input_seed_rain);
	
	S.resetState();
    S.initialize();

    S.print();

	vector <double> times = generateDefaultCohortSchedule(105.32);
	for (auto t : times) cout << t << " "; cout << endl;

	
	SolverIO<PlantModel, LightEnvironment> sio;
	sio.S = &S;
	sio.openStreams();

	ofstream fli("light_profile_ind_plant.txt");
	
	vector <vector<double>> seeds_out(S.n_species());

	for (size_t i=0; i < times.size(); ++i){

		S.step_to(times[i]);		
		
		vector<double> seeds = S.newborns_out();
		for (int s=0; s< S.n_species(); ++s){
			double S_D = 0.25;
			seeds_out[s].push_back(seeds[s] * S_D * env.patch_age_density(times[i]));
		}

		cout << times[i] << " " << S.get_species(0)->xsize() << " " << seeds[0] << " " << env.light_profile.npoints << " | " << M.nrc << " " << M.ndc << " " << M.nbc <<"\n";

		vector<double> xl = seq(0, 20, 200);
		for (auto h : xl) fli << env.canopy_openness(h) << "\t";
		fli << endl;

		sio.writeState();

	}
	
	fli.close();
	sio.closeStreams();
	cout << "derivative computations in g/m/f functions: " << M.nrc << " " << M.ndc << " " << M.nbc << endl;

	for (int s=0; s< S.n_species(); ++s){
		//auto spp = S.get_species(s);
		//auto iset = spp->get_iterators(S.state);
		//auto& itf = iset.get("fec");
		//vector <double> fec_vec;
		//fec_vec.reserve(spp->xsize());
		//iset.rbegin();
		//for (int i=0; !iset.rend(); --iset, ++i){
		//    double patch_age_density = env.patch_age_density(times[i]);
		//    double S_D = 0.25;
		//    double output_seeds = spp->mod->input_seed_rain * S_D * patch_age_density * (*itf);
		//    //cout << times[i] << " " << M.input_seed_rain << " " << S_D << " " << patch_age_density << " " << (*itf) << " | " << output_seeds << endl;
		//    fec_vec.push_back(output_seeds);
		//}
		//cout << "Seed rain for Species " << s << " = " << pn::integrate_trapezium(times, fec_vec) << endl;
		cout << "Seed rain for Species " << s << " (new method) = " << pn::integrate_trapezium(times, seeds_out[s]) << endl;

	}

}

