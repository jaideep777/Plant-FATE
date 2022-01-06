#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include <solver.h>
#include "pspm_interface.h"
#include "trait_reader.h"

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

	void openStreams(vector<string> varnames, string dir = "."){
		varnames.insert(varnames.begin(), "u");
		varnames.insert(varnames.begin(), "X");
		
		for (int s=0; s < S->species_vec.size(); ++s){
			auto spp = S->species_vec[s];
			vector<ofstream> spp_streams;
			
			for (int i=0; i<varnames.size(); ++i){
				stringstream sout;
				sout << dir << "/species_" << s << "_" << varnames[i] << ".txt";
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
				int is = 0;
				streams[s][is++] << C.x << "\t";
				streams[s][is++] << C.u << "\t";
				streams[s][is++] << C.geometry.height << "\t";
				streams[s][is++] << C.geometry.lai << "\t";
				streams[s][is++] << C.rates.dmort_dt << "\t";
				streams[s][is++] << C.state.seed_pool << "\t";
				streams[s][is++] << C.rates.rgr << "\t";
				streams[s][is++] << C.res.gpp/C.geometry.crown_area << "\t";

			}
			
			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << endl; //"\n";
		}
	}
};


int main(){


	PSPM_Dynamic_Environment E;
	E.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	E.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	E.init();
	E.print(0);
	E.use_ppa = true;

	io::Initializer I("tests/params/p.ini");
	I.readFile();
	string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
	string command = "mkdir -p " + out_dir;
	system(command.c_str());

	Solver S(SOLVER_IFMU, "rk45ck");
    S.control.ode_ifmu_stepsize = 0.0833333;
	S.control.ifmu_centered_grids = false; //true;
    S.use_log_densities = true;
	S.setEnvironment(&E);
	
	TraitsReader Tr;
	Tr.readFromFile("tests/data/early_late.csv");
	Tr.print();

	for (int i=0; i<2; ++i){
		PSPM_Plant p1;
		p1.initParamsFromFile("tests/params/p.ini");
		p1.traits.species_name = Tr.species[i].species_name;
		p1.traits.lma = Tr.species[i].lma;
		p1.traits.wood_density = Tr.species[i].wood_density;
		p1.traits.hmat = Tr.species[i].hmat;
		p1.coordinateTraits();

		((plant::Plant*)&p1)->print();
		
		//p1.geometry.set_lai(p1.par.lai0);	
		p1.set_size(0.01);
		
		Species<PSPM_Plant>* spp = new Species<PSPM_Plant>(p1);

		S.addSpecies(30, 0.01, 10, true, spp, 3, -1);
		
		//	S.addSpecies(vector<double>(1, p1.geometry.get_size()), &spp, 3, 1);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	}
	S.resetState();
	S.initialize();

	for (auto spp : S.species_vec) spp->setU(0, 1);
	S.copyCohortsToState();

	S.print();
	S.current_time = 0;
//	S.control.update_cohorts = false;


	SolverIO sio;
	sio.S = &S;
	sio.openStreams({"height", "lai", "mort", "seeds", "g", "gpp"}, out_dir);


	
//	S.step_to(0.1);
//	S.print();
//	for (auto y : S.state) cout << y << "\t"; cout << "\n";

//	ofstream fout("fmu_PlantFATE.txt");
	ofstream fzst(string(out_dir + "/z_star.txt").c_str());
	ofstream fco(string(out_dir + "/canopy_openness.txt").c_str());
	ofstream fseed(string(out_dir + "/seeds.txt").c_str());
	ofstream fabase(string(out_dir + "/basal_area.txt").c_str());
	ofstream flai(string(out_dir + "/LAI.txt").c_str());
	double t_clear = 20000;
	for (double t=0.1; t <= 800; t=t+2) {
		cout << "t = " << t << endl; //"\t";
		S.step_to(t);
		
		vector<double> seeds = S.newborns_out();
//		for (int s=0; s< S.species_vec.size(); ++s){
//			double S_D = 0.25;
//			seeds_out[s].push_back(seeds[s] * env.patch_age_density(t));
//		}
		fseed << t << "\t";
		for (int i=0; i<S.n_species(); ++i) fseed << seeds[i] << "\t";
		fseed << "\n";
		
		vector<double> basal_area(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			basal_area[k] = S.integrate_wudx_above([&S,k](int i, double t){
										  double D = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry.diameter;
										  return M_PI*D*D/4;
										}, t, 0.1, k);
		
		fabase << t << "\t";
		for (int i=0; i<S.n_species(); ++i) fabase << basal_area[i] << "\t";
		fabase << "\n";
		
		double comm_lai = 0;
		for (int k=0; k<S.n_species(); ++k)
			comm_lai += S.integrate_x([&S,k](int i, double t){
										  auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry;
										  return p.crown_area*p.lai;
										}, t, k);
		
		flai << t << "\t" << comm_lai << "\n";

//		cout << t << " " << S.species_vec[0]->xsize() << " ";
//		for (int i=0; i<S.n_species(); ++i) cout << seeds[i] << " ";
//		cout << " | n_lp = " << env.light_profile.npoints << " | PPA: nl = " << env.n_layers << ", z* = " << env.z_star[0] << "\n";

//		vector<double> xl = myseq(0, 20, 1000);
//		for (auto h : xl) fli << env.LA_above_z(t,h,&S) << "\t";
//		fli << endl;
//		
		fzst << t << "\t";
		for (auto z : E.z_star) fzst << z << "\t";
		fzst << endl;
		
		fco << t << "\t";
		for (auto z : E.canopy_openness) fco << z << "\t";
		fco << endl;
			
		sio.writeState();
	
		// clear patch after 50 year	
		if (t >= t_clear){
			for (auto spp : S.species_vec) 
				for (int i=1; i<spp->xsize(); ++i) spp->setU(i, 0);
			S.copyCohortsToState();
			t_clear = t + 200 + 100*(2*double(rand())/RAND_MAX - 1);
		}
	}
	fco.close();
	fseed.close();
	fzst.close();
	fabase.close();
	flai.close();

	// free memory associated
	for (auto s : S.species_vec) delete static_cast<Species<PSPM_Plant>*>(s);
}

