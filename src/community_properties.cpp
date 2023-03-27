#include "community_properties.h"

using namespace std;


void SpeciesProps::resize(int n){
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

SpeciesProps& SpeciesProps::operator /= (double s){
	
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

SpeciesProps& SpeciesProps::operator += (const SpeciesProps &s){
	
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

bool SpeciesProps::isResident(Species_Base * spp){
	return static_cast<MySpecies<PSPM_Plant>*>(spp)->isResident;
}

void SpeciesProps::update(double t, Solver &S){
		
	n_ind_vec.clear();
	n_ind_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		n_ind_vec[k] = S.integrate_x([&S,k](int i, double t){
											return 1;
									}, t, k);
	n_ind = std::accumulate(n_ind_vec.begin(), n_ind_vec.end(), 0.0);

	biomass_vec.clear();
	biomass_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		biomass_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											return p.get_biomass();
									}, t, 0.1, k);
	biomass = std::accumulate(biomass_vec.begin(), biomass_vec.end(), 0.0);

	ba_vec.clear();
	ba_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		ba_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											double D = p.geometry.diameter_at_height(1.3, p.traits);
											return M_PI*D*D/4;
									}, t, 0.1, k);
	ba = std::accumulate(ba_vec.begin(), ba_vec.end(), 0.0);

	canopy_area_vec.clear();
	canopy_area_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		canopy_area_vec[k] = S.integrate_x([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry;
											return p.crown_area;
									}, t, k);
	canopy_area = std::accumulate(canopy_area_vec.begin(), canopy_area_vec.end(), 0.0);
	

	height_vec.clear();
	height_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		height_vec[k] = S.integrate_x([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											return p.geometry.height;
									}, t, k);
									
	for (int k=0; k<S.n_species(); ++k) height_vec[k] /= n_ind_vec[k];


	hmat = 0;
	// for (int k=0; k<S.n_species(); ++k)
	// 	hmat += S.integrate_x([&S,k](int i, double t){
	// 								      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
	// 								      return p.traits.hmat;
	// 								}, t, k);
	// hmat /= n_ind;
	hmat_vec.clear();
	hmat_vec.resize(S.n_species(), 0);
	// for (int k=0; k<S.n_species(); ++k) hmat_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.hmat;


	lma = 0;
	// for (int k=0; k<S.n_species(); ++k)
	// 	lma += S.integrate_x([&S,k](int i, double t){
	// 								      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
	// 								      return p.traits.lma;
	// 								}, t, k);
	// lma /= n_ind;
	lma_vec.clear();
	lma_vec.resize(S.n_species(), 0);
	// for (int k=0; k<S.n_species(); ++k) lma_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.lma;

	wd = 0;
	// for (int k=0; k<S.n_species(); ++k)
	// 	wd += S.integrate_x([&S,k](int i, double t){
	// 								      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
	// 								      return p.traits.wood_density;
	// 								}, t, k);
	// wd /= n_ind;
	wd_vec.clear();
	wd_vec.resize(S.n_species(), 0);
	// for (int k=0; k<S.n_species(); ++k) wd_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.wood_density;

	p50 = 0;
	// for (int k=0; k<S.n_species(); ++k)
	// 	p50 += S.integrate_x([&S,k](int i, double t){
	// 								      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
	// 								      return p.traits.p50_xylem;
	// 								}, t, k);
	// p50 /= n_ind;
	p50_vec.clear();
	p50_vec.resize(S.n_species(), 0);
	// for (int k=0; k<S.n_species(); ++k) p50_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.p50_xylem;

	// This should not be used
	gs = 0;
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		gs += S.integrate_x([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											return p.res.gs_avg * p.geometry.crown_area;
									}, t, k);
	gs /= canopy_area;

	vcmax_vec.clear();
	vcmax_vec.resize(S.n_species(), 0);
	for (int k=0; k<S.n_species(); ++k)
		if (isResident(S.species_vec[k]))
		vcmax_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
											auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											return p.res.vcmax_avg * p.geometry.crown_area;
									}, t, 0.1, k);
	vcmax = std::accumulate(vcmax_vec.begin(), vcmax_vec.end(), 0.0);
	vcmax /= canopy_area;
}


SpeciesProps operator + (SpeciesProps lhs, SpeciesProps &rhs){
	lhs += rhs;
	return lhs;
}


EmergentProps& EmergentProps::operator /= (double s){
	
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
	
	transform(lai_vert.begin(), lai_vert.end(), lai_vert.begin(), [s](const double &c){ return c/s; });
	
	return *this;
}

EmergentProps & EmergentProps::operator += (const EmergentProps &s){
	
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

	transform(lai_vert.begin(), lai_vert.end(), s.lai_vert.begin(), lai_vert.begin(), std::plus<double>());

	return *this;
}

bool EmergentProps::isResident(Species_Base * spp){
	return static_cast<MySpecies<PSPM_Plant>*>(spp)->isResident;
}


void EmergentProps::update(double t, Solver &S){
	gpp = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.gpp;});
	npp = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.npp;});
	trans = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.trans;});
	resp_auto = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.rleaf + p->res.rroot + p->res.rstem;});
	lai = integrate_prop(t, S, [](const PSPM_Plant* p){return p->geometry.crown_area*p->geometry.lai;});
	leaf_mass = integrate_prop(t, S, [](const PSPM_Plant* p){return p->geometry.leaf_mass(p->traits);});
	stem_mass = integrate_prop(t, S, [](const PSPM_Plant* p){return p->geometry.stem_mass(p->traits);});
	croot_mass = integrate_prop(t, S, [](const PSPM_Plant* p){return p->geometry.coarse_root_mass(p->traits);});
	froot_mass = integrate_prop(t, S, [](const PSPM_Plant* p){return p->geometry.root_mass(p->traits);});
	gs = (trans*55.55/365/86400)/1.6/(static_cast<PSPM_Dynamic_Environment*>(S.env)->clim.vpd/1.0325e5);
	//     ^ convert kg/m2/yr --> mol/m2/s

	double tleaf_comm = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.tleaf;});
	double troot_comm = integrate_prop(t, S, [](const PSPM_Plant* p){return p->res.troot;});
	cc_est = (tleaf_comm + troot_comm + resp_auto)/tleaf_comm;

	// LAI vertical profile
	lai_vert.clear();
	lai_vert.resize(25, 0);
	for (int iz=0; iz<25; ++iz)
		for (int k=0; k<S.n_species(); ++k)
			if (isResident(S.species_vec[k]))
			lai_vert[iz] += S.integrate_x([&S,k,iz](int i, double t){
									auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
									return p.geometry.crown_area_above(iz,p.traits)*p.geometry.lai;
							}, t, k);

}



EmergentProps operator + (EmergentProps lhs, EmergentProps &rhs){
	lhs += rhs;
	return lhs;
}



void SolverIO::openStreams(std::string dir, io::Initializer &I){

	cohort_props_out.open(dir + "/cohort_props.txt");
	cohort_props_out << "t\tspeciesID\tcohortID\t";
	for (auto vname : varnames) cohort_props_out << vname << "\t";
	cohort_props_out << std::endl;

	size_dists_out.open(dir + "/size_distributions.txt");

	// varnames.insert(varnames.begin(), "u");
	// varnames.insert(varnames.begin(), "X");
	
	// for (int s=0; s < S->species_vec.size(); ++s){
	// 	auto spp = S->species_vec[s];
	// 	std::vector<std::ofstream> spp_streams;
		
	// 	for (int i=0; i<varnames.size(); ++i){
	// 		std::stringstream sout;
	// 		sout << dir << "/species_" << s << "_" << varnames[i] << ".txt";
	// 		cout << sout.str() << endl;
	// 		std::ofstream fout(sout.str().c_str());
	// 		assert(fout);
	// 		spp_streams.push_back(std::move(fout));
	// 	}
	// 	streams.push_back(std::move(spp_streams));
	// }

	fzst.open(std::string(dir + "/z_star.txt").c_str());
	fco.open(std::string(dir + "/canopy_openness.txt").c_str());
	// fseed.open(std::string(dir + "/seeds.txt").c_str());
	// fabase.open(std::string(dir + "/basal_area.txt").c_str());
	flai.open(std::string(dir + "/lai_profile.txt").c_str());
	foutd.open(std::string(dir + "/" + I.get<std::string>("emgProps")).c_str());
	fouty.open(std::string(dir + "/" + I.get<std::string>("cwmAvg")).c_str());
	fouty_spp.open(std::string(dir + "/" + I.get<std::string>("cwmperSpecies")).c_str());
	ftraits.open(std::string(dir + "/" + I.get<std::string>("traits")).c_str());

	foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\tVCMAX\tCCEST\n";
	fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\tSEEDS\n";
	ftraits << "YEAR\tSPP\tRES\tLMA\tWD\tHMAT\tP50X\tr0_last\tr0_avg\tr0_exp\tr0_cesaro\n";

}

void SolverIO::closeStreams(){
	// for (int s=0; s<streams.size(); ++s){
	// 	for (int j=0; j<streams[s].size(); ++j){
	// 		streams[s][j].close();
	// 	}
	// }
	cohort_props_out.close();
	size_dists_out.close();
	
	fzst.close();
	fco.close();
	// fseed.close();
	// fabase.close();
	flai.close();
	foutd.close();
	fouty.close();
	fouty_spp.close();
	ftraits.close();

}

void SolverIO::writeState(double t, SpeciesProps& cwm, EmergentProps& props){
	for (int s=0; s < S->species_vec.size(); ++s){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[s]);

		// for (int i=0; i<streams[s].size(); ++i) streams[s][i] << t << "\t";

		std::vector<double> breaks = my_log_seq(0.01, 10, 100);
		std::vector<double> dist = S->getDensitySpecies(s, breaks);
		//cout << "here: " << breaks.size() << " " << dist.size() << endl;

		if (spp->isResident){
			size_dists_out << t << "\t" << spp->species_name << "\t";
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
								<< spp->species_name << "\t"  // use name instead of index s becuase it is unique and order-insensitive
								<< j << "\t"
								<< C.geometry.height << "\t"
								<< C.geometry.lai << "\t"
								<< C.rates.dmort_dt << "\t"
								<< C.rates.dseeds_dt << "\t"
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
			<< cwm.vcmax << "\t"
			<< props.cc_est << std::endl;
	
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
			<< cwm.p50  << std::endl;
	
	for (int k=0; k<S->species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[k]);
		fouty_spp 
				<< int(t) << "\t"
				<< spp->species_name  << "\t" // use name instead of index k becuase it is unique and order-insensitive
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
				<< spp->seeds_hist.get()  << "\n";
	}

	for (int k=0; k<S->species_vec.size(); ++k){
		auto spp = static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[k]);
		ftraits 
				<< t << "\t"
				<< spp->species_name  << "\t" // use name instead of index k becuase it is unique and order-insensitive
				<< spp->isResident << "\t";
		std::vector<double> v = spp->get_traits();
		for (auto vv : v)
		ftraits << vv << "\t";
		ftraits << spp->getCohort(-1).traits.hmat << "\t"
		        << spp->getCohort(-1).traits.p50_xylem << "\t";
		ftraits << spp->r0_hist.get_last() << "\t"
				<< spp->r0_hist.get() << "\t"
				<< spp->r0_hist.get_exp(0.02) << "\t"
				<< 0 << "\n"; //spp->r0_hist.get_cesaro() << "\n";
	}

	ftraits.flush();
	fouty_spp.flush();

	flai << t << "\t";
	for (int i=0; i<props.lai_vert.size(); ++i) flai << props.lai_vert[i] << "\t";
	flai << std::endl;

	// // FIXME: Delete this
	// fseed << t << "\t";
	// for (int i=0; i<S->n_species(); ++i) fseed << static_cast<MySpecies<PSPM_Plant>*>(S->species_vec[i])->seeds_hist.get() << "\t";
	// fseed << endl;
	
	// fabase << t << "\t";
	// for (int i=0; i<S->n_species(); ++i) fabase << cwm.ba_vec[i] << "\t";
	// fabase << endl;
	

	fzst << t << "\t";
	for (auto z : static_cast<PSPM_Dynamic_Environment*>(S->env)->z_star) fzst << z << "\t";
	fzst << std::endl;
	
	fco << t << "\t";
	for (auto z : static_cast<PSPM_Dynamic_Environment*>(S->env)->canopy_openness) fco << z << "\t";
	fco << std::endl;

}


