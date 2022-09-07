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
	vector<string> varnames = {"height", "lai", "mort", "seeds", "rgr", "gpp"};

	// vector <vector<ofstream>> streams;
	ofstream cohort_props_out;
	ofstream size_dists_out;

	ofstream fzst;
	ofstream fco;
	ofstream fseed;
	ofstream fabase;
	ofstream flai;
	ofstream foutd;
	ofstream fouty;
	ofstream fouty_spp;
	ofstream ftraits;

	void openStreams(string dir, io::Initializer &I){

		cohort_props_out.open(dir + "/cohort_props.txt");
		cohort_props_out << "t\tspeciesID\tcohortID\t";
		for (auto vname : varnames) cohort_props_out << vname << "\t";
		cohort_props_out << endl;

		size_dists_out.open(dir + "/size_distributions.txt");

		// varnames.insert(varnames.begin(), "u");
		// varnames.insert(varnames.begin(), "X");
		
		// for (int s=0; s < S->species_vec.size(); ++s){
		// 	auto spp = S->species_vec[s];
		// 	vector<ofstream> spp_streams;
			
		// 	for (int i=0; i<varnames.size(); ++i){
		// 		stringstream sout;
		// 		sout << dir << "/species_" << s << "_" << varnames[i] << ".txt";
		// 		cout << sout.str() << endl;
		// 		ofstream fout(sout.str().c_str());
		// 		assert(fout);
		// 		spp_streams.push_back(std::move(fout));
		// 	}
		// 	streams.push_back(std::move(spp_streams));
		// }

		fzst.open(string(dir + "/z_star.txt").c_str());
		fco.open(string(dir + "/canopy_openness.txt").c_str());
		fseed.open(string(dir + "/seeds.txt").c_str());
		fabase.open(string(dir + "/basal_area.txt").c_str());
		flai.open(string(dir + "/lai_profile.txt").c_str());
		foutd.open(string(dir + "/" + I.get<string>("emgProps")).c_str());
		fouty.open(string(dir + "/" + I.get<string>("cwmAvg")).c_str());
		fouty_spp.open(string(dir + "/" + I.get<string>("cwmperSpecies")).c_str());
		ftraits.open(string(dir + "/" + I.get<string>("traits")).c_str());

		foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\tVCMAX\n";
		fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
		fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\tSEEDS\n";
		ftraits << "YEAR\tSPP\tRES\tLMA\tWD\tr0_last\tr0_avg\tr0_exp\tr0_cesaro\n";

	}

	void closeStreams(){
		// for (int s=0; s<streams.size(); ++s){
		// 	for (int j=0; j<streams[s].size(); ++j){
		// 		streams[s][j].close();
		// 	}
		// }
		cohort_props_out.close();
		size_dists_out.close();
		
		fzst.close();
		fco.close();
		fseed.close();
		fabase.close();
		flai.close();
		foutd.close();
		fouty.close();
		fouty_spp.close();
		ftraits.close();

	}

	void writeState(double t, vector<MovingAverager>& seeds_hist, SpeciesProps& cwm, EmergentProps& props){
		for (int s=0; s < S->species_vec.size(); ++s){
			auto spp = (MySpecies<PSPM_Plant>*)S->species_vec[s];

			// for (int i=0; i<streams[s].size(); ++i) streams[s][i] << t << "\t";

			vector<double> breaks = my_log_seq(0.01, 10, 100);
			vector<double> dist = S->getDensitySpecies(s, breaks);
			//cout << "here: " << breaks.size() << " " << dist.size() << endl;

			if (spp->isResident){
				size_dists_out << t << "\t" << s << "\t";
				for (int i=0; i<100; ++i){
					// streams[s][0] << breaks[i] << "\t";
					// streams[s][1] << dist[i] << "\t";
					size_dists_out << dist[i] << "\t";
				}
				size_dists_out << "\n";
			}

			// for (int j=0; j<spp->xsize(); ++j){
			// 	auto& C = spp->getCohort(j);
			// 	int is = 2;
			// 	//streams[s][is++] << C.x << "\t";
			// 	//streams[s][is++] << C.u << "\t";
			// 	streams[s][is++] << C.geometry.height << "\t";
			// 	streams[s][is++] << C.geometry.lai << "\t";
			// 	streams[s][is++] << C.rates.dmort_dt << "\t";
			// 	streams[s][is++] << C.state.seed_pool << "\t";
			// 	streams[s][is++] << C.rates.rgr << "\t";
			// 	streams[s][is++] << C.res.gpp/C.geometry.crown_area << "\t";

			// }
			
			// for (int i=0; i<streams[s].size(); ++i) streams[s][i] << endl; //"\n";
		
			if (spp->isResident){
				for (int j=0; j<spp->xsize()-1; ++j){
					auto& C = spp->getCohort(j);
					cohort_props_out << t << "\t" 
									<< s << "\t"
									<< j << "\t"
									<< C.geometry.height << "\t"
									<< C.geometry.lai << "\t"
									<< C.rates.dmort_dt << "\t"
									<< C.state.seed_pool << "\t"
									<< C.rates.rgr << "\t"
									<< C.res.gpp/C.geometry.crown_area << "\t";
					cohort_props_out << "\n";
				}
			}
		}

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
		
		for (int k=0; k<S->species_vec.size(); ++k){
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
				  << cwm.p50_vec[k]  << "\t"
				  << seeds_hist[k].get()  << "\n";
		}

		for (int k=0; k<S->species_vec.size(); ++k){
			auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[k]);
			ftraits 
			      << t << "\t"
				  << k << "\t"
				  << spp->isResident << "\t";
			vector<double> v = spp->get_traits();
			for (auto vv : v)
			ftraits << vv << "\t";
			ftraits << spp->r0_hist.get_last() << "\t"
				    << spp->r0_hist.get() << "\t"
				    << spp->r0_hist.get_exp(0.02) << "\t"
				    << 0 << "\n"; //spp->r0_hist.get_cesaro() << "\n";
		}

		ftraits.flush();
		fouty_spp.flush();

		flai << t << "\t";
		for (int i=0; i<props.lai_vert.size(); ++i) flai << props.lai_vert[i] << "\t";
		flai << endl;

		fseed << t << "\t";
		for (int i=0; i<S->n_species(); ++i) fseed << seeds_hist[i].get() << "\t";
		fseed << endl;
		
		fabase << t << "\t";
		for (int i=0; i<S->n_species(); ++i) fabase << cwm.ba_vec[i] << "\t";
		fabase << endl;
		
	
		fzst << t << "\t";
		for (auto z : static_cast<PSPM_Dynamic_Environment*>(S->env)->z_star) fzst << z << "\t";
		fzst << endl;
		
		fco << t << "\t";
		for (auto z : static_cast<PSPM_Dynamic_Environment*>(S->env)->canopy_openness) fco << z << "\t";
		fco << endl;

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
	E.update_co2 = false;

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

	for (auto spp : S.species_vec){
		for (auto m : static_cast<MySpecies<PSPM_Plant>*>(spp)->probes) 
			S.addSpecies(res, 0.01, 10, true, m, 3, 1e-3);
	}

	S.resetState(I.getScalar("year0"));
	S.initialize();

	S.print();

	SolverIO sio;
	sio.S = &S;
	sio.openStreams(out_dir, I);


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
			
			auto spp = static_cast<MySpecies<PSPM_Plant>*>(S.species_vec[k]);
			double r0 = seeds_hist[k].get()/S.species_vec[k]->birth_flux_in;
			// cout << "   " << k << ": " <<  S.species_vec[k]->birth_flux_in << " --> " << seeds[k] << "/" << seeds_hist[k].get() << ", r0 = " << setprecision(8) << r0 << "\n";
			
			spp->set_inputBirthFlux(seeds_hist[k].get());
			spp->r0_hist.push(t, r0);
			// spp->r0_hist.print_summary();
		}
	};

//	ofstream fout("fmu_PlantFATE.txt");

	SpeciesProps cwm;
	EmergentProps props; 

	double t_clear = 105000;
	// t is years since 2000-01-01
	double y0, yf;
	y0 = I.getScalar("year0");
	yf = I.getScalar("yearf");
	double delta_T = I.getScalar("delta_T");
	for (double t=y0; t <= yf; t=t+delta_T) {
		cout << "stepping = " << S.current_time << "-->" << t << "\t(";
		for (auto spp : S.species_vec) cout << spp->xsize() << ", ";
		cout << ")" << endl;

		S.step_to(t, after_step);
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

