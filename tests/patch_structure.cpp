#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <string>
using namespace std;

#include <solver.h>
#include "pspm_interface.h"
#include "trait_reader.h"

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
	if (tvec.size() >=1)     // JAI: added to avoid overflow warning
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

class SpeciesProps{
public:
	double n_ind=0;
	double biomass=0;
	double ba=0;
	double canopy_area=0;
	double height=0;
	double lma=0;
	double p50=0;
	double hmat=0;
	double wd=0;
	double gs=0;
	double vcmax=0;
	
	vector<double> n_ind_vec;
	vector<double> biomass_vec;
	vector<double> ba_vec;
	vector<double> canopy_area_vec;
	vector<double> height_vec;
	vector<double> vcmax_vec;
	
	vector<double> lma_vec;
	vector<double> p50_vec;
	vector<double> hmat_vec;
	vector<double> wd_vec;
	
	void resize(int n){
		n_ind_vec.resize(n);
		ba_vec.resize(n);
		biomass_vec.resize(n);
		canopy_area_vec.resize(n);
		height_vec.resize(n);
		vcmax_vec.resize(n);
		hmat_vec.resize(n);
		lma_vec.resize(n);
		wd_vec.resize(n);
		p50_vec.resize(n);
	}
	
	SpeciesProps & operator /= (double s){
		
		n_ind/=s;
		biomass/=s;
		ba/=s;
		canopy_area/=s;
		height/=s;
		vcmax/=s;
		lma/=s;
		p50/=s;
		hmat/=s;
		wd/=s;
		gs/=s;
		transform(n_ind_vec.begin(), n_ind_vec.end(), n_ind_vec.begin(), [s](const double &c){ return c/s; });
		transform(biomass_vec.begin(), biomass_vec.end(), biomass_vec.begin(), [s](const double &c){ return c/s; });
		transform(ba_vec.begin(), ba_vec.end(), ba_vec.begin(), [s](const double &c){ return c/s; });
		transform(canopy_area_vec.begin(), canopy_area_vec.end(), canopy_area_vec.begin(), [s](const double &c){ return c/s; });
		transform(height_vec.begin(), height_vec.end(), height_vec.begin(), [s](const double &c){ return c/s; });
		transform(vcmax_vec.begin(), vcmax_vec.end(), vcmax_vec.begin(), [s](const double &c){ return c/s; });
		transform(lma_vec.begin(), lma_vec.end(), lma_vec.begin(), [s](const double &c){ return c/s; });
		transform(p50_vec.begin(), p50_vec.end(), p50_vec.begin(), [s](const double &c){ return c/s; });
		transform(hmat_vec.begin(), hmat_vec.end(), hmat_vec.begin(), [s](const double &c){ return c/s; });
		transform(wd_vec.begin(), wd_vec.end(), wd_vec.begin(), [s](const double &c){ return c/s; });
		
		return *this;
	}
	
	SpeciesProps & operator += (const SpeciesProps &s){
		
		n_ind+=s.n_ind;
		biomass+=s.biomass;
		ba+=s.ba;
		canopy_area+=s.canopy_area;
		height+=s.height;
		vcmax+=s.vcmax;
		lma+=s.lma;
		p50+=s.p50;
		hmat+=s.hmat;
		wd+=s.wd;
		gs+=s.gs;
		transform(n_ind_vec.begin(), n_ind_vec.end(), s.n_ind_vec.begin(), n_ind_vec.begin(), std::plus<double>());
		transform(biomass_vec.begin(), biomass_vec.end(), s.biomass_vec.begin(), biomass_vec.begin(),std::plus<double>());
		transform(ba_vec.begin(), ba_vec.end(), s.ba_vec.begin(), ba_vec.begin(),std::plus<double>());
		transform(canopy_area_vec.begin(), canopy_area_vec.end(), s.canopy_area_vec.begin(), canopy_area_vec.begin(),std::plus<double>());
		transform(height_vec.begin(), height_vec.end(), s.height_vec.begin(), height_vec.begin(),std::plus<double>());
		transform(vcmax_vec.begin(), vcmax_vec.end(), s.vcmax_vec.begin(), vcmax_vec.begin(), std::plus<double>());
		transform(lma_vec.begin(), lma_vec.end(), s.lma_vec.begin(), lma_vec.begin(),std::plus<double>());
		transform(p50_vec.begin(), p50_vec.end(), s.p50_vec.begin(), p50_vec.begin(),std::plus<double>());
		transform(hmat_vec.begin(), hmat_vec.end(), s.hmat_vec.begin(), hmat_vec.begin(),std::plus<double>());
		transform(wd_vec.begin(), wd_vec.end(), s.wd_vec.begin(), wd_vec.begin(),std::plus<double>());
		
		return *this;
	}
	
	template<class Func>
	void integrate_prop(vector<double> &vx, double t, Solver &S, const Func &f){
		vx.clear();
		vx.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			vx[k] = S.integrate_x([&S,k,f](int i, double t){
				auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				const PSPM_Plant * pp = &p;
				return f(pp);
			}, t, k);
		
	}

	template<class Func>
	void integrate_prop_above(vector<double> &vx, double t, Solver &S, const Func &f, double xlow){
		vx.clear();
		vx.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			vx[k] = S.integrate_wudx_above([&S,k,f](int i, double t){
				auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				const PSPM_Plant * pp = &p;
				return f(pp);
			}, t, xlow, k);
		
	}
	
	template<class Func>
	void integrate_prop(double &x, double t, Solver &S, const Func &f){
		x = 0;
		for (int k=0; k<S.n_species(); ++k)
			x += S.integrate_x([&S,k,f](int i, double t){
				auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				const PSPM_Plant * pp = &p;
				return f(pp);
			}, t, k);
		
	}

	template<class Func>
	void calc_cwm_trait(double &x, vector<double> &vx, double _n_ind, double t, Solver &S, const Func &f){
		integrate_prop(x, t, S, [&f](const PSPM_Plant* p){return f(p);});
		x /= _n_ind;

		vx.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k){
			auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1);
			const PSPM_Plant * pp = &p;
			vx[k] = f(pp);
		}
	}

	void update(double t, Solver &S){
		integrate_prop(n_ind_vec, t, S, [](const PSPM_Plant* p){return 1;});
		n_ind = std::accumulate(n_ind_vec.begin(), n_ind_vec.end(), 0.0);
		
		integrate_prop(biomass_vec, t, S, [](const PSPM_Plant* p){return p->get_biomass();});
		biomass = std::accumulate(biomass_vec.begin(), biomass_vec.end(), 0.0);

		integrate_prop_above(ba_vec, t, S, [](const PSPM_Plant* p){
			double D = p->geometry.diameter; 
			return M_PI*D*D/4;
			}, 0.1);
		ba = std::accumulate(ba_vec.begin(), ba_vec.end(), 0.0);
		
		integrate_prop(canopy_area_vec, t, S, [](const PSPM_Plant* p){return p->geometry.crown_area;});
		canopy_area = std::accumulate(canopy_area_vec.begin(), canopy_area_vec.end(), 0.0);
		
		integrate_prop(height_vec, t, S, [](const PSPM_Plant* p){return p->geometry.height;}); // sum(N*h)[k]
		height = std::accumulate(height_vec.begin(), height_vec.end(), 0.0);                    // sum(N*h)
		for (int k=0; k<S.n_species(); ++k) height_vec[k] /= n_ind_vec[k];                      // sum(N*h)[k]/sum(N)[k]
		height /= n_ind;                                                                        // sum(N*h) / sum(N)
		
		calc_cwm_trait(vcmax, vcmax_vec, n_ind, t, S, [](const PSPM_Plant* p){return p->res.vcmax_avg;});

		calc_cwm_trait(hmat, hmat_vec, n_ind, t, S, [](const PSPM_Plant* p){return p->traits.hmat;});
		calc_cwm_trait(lma, lma_vec, n_ind, t, S, [](const PSPM_Plant* p){return p->traits.lma;});
		calc_cwm_trait(wd, wd_vec, n_ind, t, S, [](const PSPM_Plant* p){return p->traits.wood_density;});
		calc_cwm_trait(p50, p50_vec, n_ind, t, S, [](const PSPM_Plant* p){return p->traits.p50_xylem;});
				
		// FIXME: gs should be calc as trans/1.6D in EmergentProps
		gs = 0;
		for (int k=0; k<S.n_species(); ++k)
			gs += S.integrate_x([&S,k](int i, double t){
				auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				return p.res.gs_avg;
			}, t, k);
		gs /= n_ind;
		
	}
};

SpeciesProps operator + (SpeciesProps lhs, SpeciesProps &rhs){
	lhs += rhs;
	return lhs;
}


class EmergentProps{
public:
	double gpp=0;
	double npp=0;
	double resp_auto=0;
	double trans=0;
	double gs=0;
	double lai=0;
	double leaf_mass=0;
	double stem_mass=0;
	double croot_mass=0;
	double froot_mass=0;
	
	
	EmergentProps & operator /= (double s){
		
		gpp/=s;
		npp/=s;
		resp_auto/=s;
		trans/=s;
		gs/=s;
		lai/=s;
		leaf_mass/=s;
		stem_mass/=s;
		croot_mass/=s;
		froot_mass/=s;
		
		return *this;
	}
	
	EmergentProps & operator += (const EmergentProps &s){
		
		gpp+=s.gpp;
		npp+=s.npp;
		resp_auto+=s.resp_auto;
		trans+=s.trans;
		gs+=s.gs;
		lai+=s.lai;
		leaf_mass+=s.leaf_mass;
		stem_mass+=s.stem_mass;
		croot_mass+=s.croot_mass;
		froot_mass+=s.froot_mass;
		
		return *this;
	}
	
	template<class Func>
	void integrate_prop(double &x, double t, Solver &S, const Func &f){
		x = 0;
		for (int k=0; k<S.n_species(); ++k)
			x += S.integrate_x([&S,k,f](int i, double t){
				auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				const PSPM_Plant * pp = &p;
				return f(pp);
			}, t, k);
		
	}
	
	void update(double t, Solver &S){
		integrate_prop(gpp, t, S, [](const PSPM_Plant* p){return p->res.gpp;});
		integrate_prop(npp, t, S, [](const PSPM_Plant* p){return p->res.npp;});
		integrate_prop(trans, t, S, [](const PSPM_Plant* p){return p->res.trans;});
		integrate_prop(resp_auto, t, S, [](const PSPM_Plant* p){return p->res.rleaf + p->res.rroot + p->res.rstem;});
		integrate_prop(lai, t, S, [](const PSPM_Plant* p){return p->geometry.crown_area*p->geometry.lai;});
		integrate_prop(leaf_mass, t, S, [](const PSPM_Plant* p){return p->geometry.leaf_mass(p->traits);});
		integrate_prop(stem_mass, t, S, [](const PSPM_Plant* p){return p->geometry.stem_mass(p->traits);});
		integrate_prop(croot_mass, t, S, [](const PSPM_Plant* p){return p->geometry.coarse_root_mass(p->traits);});
		integrate_prop(froot_mass, t, S, [](const PSPM_Plant* p){return p->geometry.root_mass(p->traits);});
		gs = (trans*55.55/365/86400)/1.6/(static_cast<PSPM_Dynamic_Environment*>(S.env)->clim.vpd/1.0325e5);
		//     ^ convert kg/m2/yr --> mol/m2/s
	}

};

EmergentProps operator + (EmergentProps lhs, EmergentProps &rhs){
	lhs += rhs;
	return lhs;
}


class Patch{
public:
	
	Solver S;
	PSPM_Dynamic_Environment E;
	SpeciesProps cwm;
	EmergentProps props;
	TraitsReader Tr;
	SolverIO sio;
	ofstream fzst;
	ofstream fco;
	ofstream fseed;
	ofstream fabase;
	ofstream foutd;
	ofstream fouty;
	ofstream fouty_spp;
	double t_ret;
	double t_clear;
	
	vector<double> patch_seeds;
	double T_seed_rain_avg;
	int npatch;
	
	Patch(PSPM_SolverType _method, string ode_method) : S(_method, ode_method) {
	}
	
	void initPatch(io::Initializer &I, int patchNo){
		
		npatch = I.getScalar("nPatches");
		t_ret = I.getScalar("T_return");
		E.metFile = I.get<string>("metFile");
		E.co2File = I.get<string>("co2File");
		
		//t_clear = fmin(-log(double(rand())/RAND_MAX)*t_ret,1000);
		t_clear =1050;
		
		string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
		E.init();
		E.print(0);
		E.use_ppa = true;
		E.update_met = true;
		E.update_co2 = true;

		S.control.ode_ifmu_stepsize = 0.0833333;
		S.control.ifmu_centered_grids = false; //true;
		S.use_log_densities = true;
		S.setEnvironment(&E);
		
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
			// FIXME: Need initial size calculation from seed mass
			p1.set_size(0.01);
			Species<PSPM_Plant>* spp = new Species<PSPM_Plant>(p1);
			S.addSpecies(30, 0.01, 10, true, spp, 3, 1e-3);
			//S.addSpecies({0.01,0.01+1e-4}, spp, 3, 1e-3);
			
			//    S.addSpecies(vector<double>(1, p1.geometry.get_size()), &spp, 3, 1);
			//S.get_species(0)->set_bfin_is_u0in(true);    // say that input_birth_flux is u0
		}
		
		S.resetState(I.getScalar("year0"));
		S.initialize();
		
		for (auto spp : S.species_vec) spp->setU(0, 1);
		S.copyCohortsToState();
		
		S.print();
		sio.S = &S;
	}
	
	void initOutputFiles(io::Initializer &I, int patchNo){
		string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
		string patchpath = string(out_dir + "/Patch"+ to_string(patchNo));
		string command = "mkdir -p " + patchpath;
		int sysresult;
		sysresult = system(command.c_str());
		
		sio.openStreams({"height", "lai", "mort", "seeds", "g", "gpp"}, out_dir+"/Patch"+ to_string(patchNo));
		
		fzst.open(string(patchpath + "/z_star.txt").c_str());
		fco.open(string(patchpath + "/canopy_openness.txt").c_str());
		fseed.open(string(patchpath + "/seeds.txt").c_str());
		fabase.open(string(patchpath + "/basal_area.txt").c_str());
		foutd.open(string(patchpath + "/" + I.get<string>("emgProps")).c_str());
		fouty.open(string(patchpath + "/" + I.get<string>("cwmAvg")).c_str());
		fouty_spp.open(string(patchpath + "/" + I.get<string>("cwmperSpecies")).c_str());
		foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\n";
		fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
		fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	}
	
	void stepPatch(double t){
		S.step_to(t);
		patch_seeds = S.newborns_out(t);
	}
	
	
	void writePatchData(double t, vector<MovingAverager> seeds_hist){
		
		
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
		<< props.gs << "\t"
		<< props.trans/365 << "\t"   // kg/m2/yr --> 1e-3 m3/m2/yr --> 1e-3*1e3 mm/yr --> 1/365 mm/day
		<< props.lai << endl;
		
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
		
		fco.flush();
		fseed.flush();
		fzst.flush();
		fabase.flush();

		sio.writeState();
		
	}
	
	
	void clearPatch(io::Initializer &I,double t){
		
		if (t >= t_clear){
			for (auto spp : S.species_vec){
				for (int i=0; i<spp->xsize(); ++i){
					auto& p = (static_cast<Species<PSPM_Plant>*>(spp))->getCohort(i);
					p.geometry.lai = p.par.lai0;
					double u_new = 0; //spp->getU(i) * double(rand())/RAND_MAX;
					spp->setU(i, u_new);
				}
			}
			S.copyCohortsToState();
			double t_int = -log(double(rand())/RAND_MAX) * t_ret;
			t_clear = t + fmin(t_int, 1000);
		}
	}
	
	void closePatch(){
		fco.close();
		fseed.close();
		fzst.close();
		fabase.close();
		foutd.close();
		fouty.close();
		for (auto s : S.species_vec) delete static_cast<Species<PSPM_Plant>*>(s);
	}
	
};


// AMZ-FACE TODO:
//	vcmax is too low (16) vs obs (40)
//	gs is too high (0.3 / 0.6) vs obs (0.15)
	

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
	
	int nspp_global = I.getScalar("nSpecies");
	int npatches = I.getScalar("nPatches");
	
	// Initialising Vector of Patch Object Pointers
	vector<Patch*> pa;
	// Initialising Moving Average function for Seed history of Ecosystem
	vector<MovingAverager> seeds_hist;
	
	// Initialising Individual Patches
	for (int i=0;i<npatches;i++){
		Patch *patch = new Patch(SOLVER_IFMU, "rk45ck");
		pa.push_back(patch);
	}
	for (int i=0; i<npatches; ++i){
		pa[i]->initPatch(I, i);
		pa[i]->initOutputFiles(I, i);
	}
	
	// Creating Output files for Ecosystem
	ofstream fseed(string(out_dir + "/seeds.txt").c_str());
	ofstream fabase(string(out_dir + "/basal_area.txt").c_str());
	ofstream foutd(string(out_dir + "/" + I.get<string>("emgProps")).c_str());
	ofstream fouty(string(out_dir + "/" + I.get<string>("cwmAvg")).c_str());
	ofstream fouty_spp(string(out_dir + "/" + I.get<string>("cwmperSpecies")).c_str());
	foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\n";
	fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\tVCM\n";
	fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	
	
	double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
	seeds_hist.resize(nspp_global);
	for (auto& M : seeds_hist) M.set_interval(T_seed_rain_avg);
	
	
	double y0, yf;
	y0 = I.getScalar("year0");
	yf = I.getScalar("yearf");
	
	// Simulation Start Point
	for (double t=y0; t <= yf; t=t+1){
		cout << "t = " << t << endl;
		SpeciesProps cwmEcosystem;
		EmergentProps propsEcosystem;
		
		cwmEcosystem.resize(nspp_global);
		
		// step all patches
		for(int i=0;i<npatches;i++){
			pa[i]->stepPatch(t);
		}
		
		// calculate average seed output from all patches 
		vector<double> total_seeds(nspp_global,0);
		for(int k=0; k<nspp_global; ++k){
			for(int i=0;i<npatches;++i){
				total_seeds[k] = total_seeds[k] + pa[i]->patch_seeds[k];
			}
		}
		for (auto& x : total_seeds) x/=npatches;
		// Push total seeds to history and update patch seed inputs
		for (int k=0; k<nspp_global; ++k){
			seeds_hist[k].push(t, total_seeds[k]);
			for(int i=0;i<npatches;i++){
				pa[i]->S.species_vec[k]->set_inputBirthFlux(seeds_hist[k].get());
			}
		}

		// output patch data
		for(int i=0;i<npatches;i++){
			pa[i]->writePatchData(t, seeds_hist);
		}
		
		// implement disturbance - clear patches if disturbance interval has been reached
		for(int i=0;i<npatches;i++){
			pa[i]->clearPatch(I,t);
		}

		// Averaging SpeciesProps and EmergentProps values for entire Ecosystem
		for(int i=0; i<npatches; i++){
			cwmEcosystem += pa[i]->cwm;
			propsEcosystem += pa[i]->props;
		}
		cwmEcosystem/=npatches;
		propsEcosystem/=npatches;
		
		// Storing the data into file streams for entire Ecosystem
		foutd << int(t) << "\t"
		<< (t-int(t))*365 << "\t"
		<< propsEcosystem.gpp*0.5/365*1000 << "\t"
		<< propsEcosystem.npp*0.5/365*1000 << "\t"
		<< propsEcosystem.resp_auto*0.5/365*1000 << "\t"  // gC/m2/d
		<< propsEcosystem.leaf_mass*1000*0.5 << "\t"
		<< propsEcosystem.stem_mass*1000*0.5 << "\t"
		<< propsEcosystem.croot_mass*1000*0.5 << "\t"
		<< propsEcosystem.froot_mass*1000*0.5 << "\t"
		<< (propsEcosystem.croot_mass+propsEcosystem.froot_mass)*1000*0.5 << "\t" // gC/m2
		<< propsEcosystem.gs << "\t"
		<< propsEcosystem.trans/365 << "\t"   // kg/m2/yr --> 1e-3 m3/m2/yr --> 1e-3*1e3 mm/yr --> 1/365 mm/day
		<< propsEcosystem.lai << endl;
		
		fouty << int(t) << "\t"
		<< -9999  << "\t"
		<< cwmEcosystem.n_ind << "\t"
		<< -9999  << "\t"
		<< cwmEcosystem.height  << "\t"
		<< cwmEcosystem.hmat  << "\t"
		<< cwmEcosystem.canopy_area  << "\t"   // m2/m2
		<< cwmEcosystem.ba  << "\t"            // m2/m2
		<< cwmEcosystem.biomass  << "\t"       // kg/m2
		<< cwmEcosystem.wd  << "\t"
		<< -9999  << "\t"
		<< 1/cwmEcosystem.lma  << "\t"
		<< cwmEcosystem.p50  << "\t"
		<< cwmEcosystem.vcmax << endl;
		
		for (int k=0; k<nspp_global; ++k){
			fouty_spp
			<< int(t) << "\t"
			<< k  << "\t"
			<< cwmEcosystem.n_ind_vec[k] << "\t"
			<< -9999  << "\t"
			<< cwmEcosystem.height_vec[k]  << "\t"
			<< cwmEcosystem.hmat_vec[k]  << "\t"
			<< cwmEcosystem.canopy_area_vec[k]  << "\t"   // m2/m2
			<< cwmEcosystem.ba_vec[k]  << "\t"            // m2/m2
			<< cwmEcosystem.biomass_vec[k]  << "\t"       // kg/m2
			<< cwmEcosystem.wd_vec[k]  << "\t"
			<< -9999  << "\t"
			<< 1/cwmEcosystem.lma_vec[k]  << "\t"
			<< cwmEcosystem.p50_vec[k]  << "\n";
		}
		
		fseed << t << "\t";
		for (int i=0; i<nspp_global; ++i) fseed << seeds_hist[i].get() << "\t";
		fseed << "\n";
		
		fabase << t << "\t";
		for (int i=0; i<nspp_global; ++i) fabase << cwmEcosystem.ba_vec[i] << "\t";
		fabase << "\n";
		
		fseed.flush();
		fabase.flush();
	}
	
	fseed.close();
	fabase.close();
	foutd.close();
	fouty.close();
	fouty_spp.close();
	for (int i=0; i<npatches; ++i){
		pa[i]->closePatch();
		delete pa[i];
	}
}



