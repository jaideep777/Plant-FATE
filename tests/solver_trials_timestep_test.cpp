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

inline double runif(double rmin=0, double rmax=1){
	double r = double(rand())/RAND_MAX; 
	return rmin + (rmax-rmin)*r;
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
				assert(fout);
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
			auto spp = (MySpecies<PSPM_Plant>*)S->species_vec[s];

			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << S->current_time << "\t";

			// vector<double> breaks = my_log_seq(0.01, 10, 100);
			// vector<double> dist = S->getDensitySpecies(s, breaks);
			// //cout << "here: " << breaks.size() << " " << dist.size() << endl;

			// for (int i=0; i<100; ++i){
			// 	streams[s][0] << breaks[i] << "\t";
			// 	streams[s][1] << dist[i] << "\t";
			// }

			for (int j=0; j<spp->xsize(); ++j){
				auto& C = spp->getCohort(j);
				int is = 2;
				//streams[s][is++] << C.x << "\t";
				//streams[s][is++] << C.u << "\t";
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

	string paramsFile = "tests/params/p.ini";
	io::Initializer I(paramsFile);
	I.readFile();
	string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
	string command = "mkdir -p " + out_dir;
	string command2 = "cp tests/params/p.ini " + out_dir;
	int sysresult;
	sysresult = system(command.c_str());
	sysresult = system(command2.c_str());

	PSPM_Dynamic_Environment E;
	E.metFile = I.get<string>("metFile");
	E.co2File = I.get<string>("co2File");
	E.init();
	E.print(0);
	E.use_ppa = true;
	E.update_met = true;
	E.update_co2 = true;

	string solver_method = I.get<string>("solver");
	Solver S(solver_method, "rk45ck");
	S.control.abm_n0 = 20;
    S.control.ode_ifmu_stepsize = I.getScalar("timestep"); //0.02; //0.0833333;
	S.control.ifmu_centered_grids = false; //true;
	S.control.ifmu_order = 1;
	S.control.ebt_ucut = 1e-7;
    S.use_log_densities = true;
	S.setEnvironment(&E);
	
	TraitsReader Tr;
	Tr.readFromFile(I.get<string>("traitsFile"));
	Tr.print();

	int nspp = I.getScalar("nSpecies");
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

		int res = I.getScalar("resolution");
		S.addSpecies(res, 0.01, 10, true, spp, 3, 1e-3);
		//S.addSpecies({0.01, 0.0100001}, spp, 3, 1e-3);
		
		//	S.addSpecies(vector<double>(1, p1.geometry.get_size()), &spp, 3, 1);
		//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	}
	S.resetState(I.getScalar("year0"));
	S.initialize();

	// for (auto spp : S.species_vec) spp->setU(0, 1);
	// S.copyCohortsToState();

	S.print();
//	S.control.update_cohorts = false;

	SolverIO sio;
	sio.S = &S;
	sio.openStreams({"height", "lai", "mort", "seeds", "rgr", "gpp"}, out_dir);



//	S.step_to(0.1);
//	S.print();
//	for (auto y : S.state) cout << y << "\t"; cout << "\n";
	
	double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
	vector<MovingAverager> seeds_hist(S.species_vec.size());
	for (auto& M : seeds_hist) M.set_interval(T_seed_rain_avg);

	auto after_step = [&S, &seeds_hist](double t){
		vector<double> seeds = S.newborns_out(t);
		for (int k=0; k<S.species_vec.size(); ++k){
			seeds_hist[k].push(t, seeds[k]);
			if (seeds[k] < 0){
				cout << "seeds[" << k << "] = " << seeds[k] << endl;
				S.print();
				seeds_hist[k].print_summary(); cout.flush();
			}
			S.species_vec[k]->set_inputBirthFlux(seeds_hist[k].get());
		}
	};

//	ofstream fout("fmu_PlantFATE.txt");
	ofstream fzst(string(out_dir + "/z_star.txt").c_str());
	assert(fzst);
	ofstream fco(string(out_dir + "/canopy_openness.txt").c_str());
	ofstream flai(string(out_dir + "/lai_profile.txt").c_str());
	ofstream fseed(string(out_dir + "/seeds.txt").c_str());
	ofstream fabase(string(out_dir + "/basal_area.txt").c_str());
//	ofstream flai(string(out_dir + "/LAI.txt").c_str());
//	ofstream fcwmt(string(out_dir + "/cwmt.txt").c_str());
	ofstream foutd(string(out_dir + "/" + I.get<string>("emgProps")).c_str());
	ofstream fouty(string(out_dir + "/" + I.get<string>("cwmAvg")).c_str());
	ofstream fouty_spp(string(out_dir + "/" + I.get<string>("cwmperSpecies")).c_str());
	
	foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\tVCMAX\n";
	fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	double t_clear = 105000;
	// t is years since 2000-01-01
	double y0, yf;
	y0 = I.getScalar("year0");
	yf = I.getScalar("yearf");
	double delta_T = I.getScalar("delta_T");
	for (double t=y0; t <= yf; t=t+delta_T) {
		cout << "t = " << t << "\t(";
		for (auto spp : S.species_vec) cout << spp->xsize() << ", ";
		cout << ")" << endl;

		S.step_to(t, after_step);
		//S.print(); cout.flush();

		SpeciesProps cwm;
		EmergentProps props; 
		
		cwm.update(t, S);
		props.update(t, S);
		
		foutd << int(t) << "\t"
			  << (t-int(t))*365 << "\t"
			  << props.gpp*0.5/365*1000 << "\t"
			  << props.npp*0.5/365*1000 << "\t"
			  << props.resp_auto*0.5/365*1000 << "\t"  // gC/m2/d
			  << props.leaf_mass*1000*0.5 << "\t"     
			  << props.stem_mass*1000*0.5 << "\t"
			  << props.croot_mass*1000*0.5 << "\t"
			  << props.froot_mass*1000*0.5 << "\t"
			  << (props.croot_mass+props.froot_mass)*1000*0.5 << "\t" // gC/m2
			  << cwm.gs << "\t"
			  << props.trans/365 << "\t"   // kg/m2/yr --> 1e-3 m3/m2/yr --> 1e-3*1e3 mm/yr --> 1/365 mm/day  
			  << props.lai << "\t"
			  << cwm.vcmax << endl;
		
		fouty << int(t) << "\t"
		      << -9999  << "\t"
		      << cwm.n_ind << "\t"
		      << -9999  << "\t"
		      << cwm.height  << "\t"
		      << cwm.hmat  << "\t"
		      << cwm.canopy_area  << "\t"   // m2/m2
		      << cwm.ba  << "\t"            // m2/m2
		      << cwm.biomass  << "\t"       // kg/m2
		      << cwm.wd  << "\t"
		      << -9999  << "\t"
		      << 1/cwm.lma  << "\t"
		      << cwm.p50  << endl;
		
		for (int k=0; k<S.species_vec.size(); ++k){
			fouty_spp 
			      << int(t) << "\t"
				  << k  << "\t"
				  << cwm.n_ind_vec[k] << "\t"
				  << -9999  << "\t"
				  << cwm.height_vec[k]  << "\t"
				  << cwm.hmat_vec[k]  << "\t"
				  << cwm.canopy_area_vec[k]  << "\t"   // m2/m2
				  << cwm.ba_vec[k]  << "\t"            // m2/m2
				  << cwm.biomass_vec[k]  << "\t"       // kg/m2
				  << cwm.wd_vec[k]  << "\t"
				  << -9999  << "\t"
				  << 1/cwm.lma_vec[k]  << "\t"
				  << cwm.p50_vec[k]  << "\n";
		}
	
		flai << t << "\t";
		for (int i=0; i<props.lai_vert.size(); ++i) flai << props.lai_vert[i] << "\t";
		flai << "\n";
		//fouty << t << -9999 << cwm.n_ind << -9999 << 
		

		fseed << t << "\t";
		for (int i=0; i<S.n_species(); ++i) fseed << seeds_hist[i].get() << "\t";
		fseed << "\n";
		
		fabase << t << "\t";
		for (int i=0; i<S.n_species(); ++i) fabase << cwm.ba_vec[i] << "\t";
		fabase << "\n";
		
	
		fzst << t << "\t";
		for (auto z : E.z_star) fzst << z << "\t";
		fzst << endl;
		
		fco << t << "\t";
		for (auto z : E.canopy_openness) fco << z << "\t";
		fco << endl;
			
		sio.writeState();
	
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
		
		fco.flush();
		fseed.flush();
		fzst.flush();
		fabase.flush();
//		flai.flush();
//		fcwmt.flush();

	}
	
	S.print();

	fco.close();
	fseed.close();
	fzst.close();
	fabase.close();
//	flai.close();
//	fcwmt.close();
	foutd.close();
	fouty.close();
	// free memory associated
	for (auto s : S.species_vec) delete static_cast<MySpecies<PSPM_Plant>*>(s); 
}

