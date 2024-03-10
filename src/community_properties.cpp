#include "community_properties.h"
#include "light_environment.h"
#include "plantfate_patch.h"

using namespace std;

namespace pfate{

void CommunityProperties::resize(int n){
	species.n_ind_vec.resize(n);
	species.biomass_vec.resize(n);
	species.basal_area_vec.resize(n);
	species.canopy_area_vec.resize(n);
	species.height_vec.resize(n);
}

bool CommunityProperties::isResident(Species_Base * spp){
	return static_cast<AdaptiveSpecies<PSPM_Plant>*>(spp)->isResident;
}

void CommunityProperties::update(double t, Patch &P){

	Solver &S = P.S;

	// multiplier to convert unit_t-1 --> day-1
	double m1 = 1/P.par0.days_per_tunit;

	// multiplier to convert kg biomass --> kg C
	double m2 = 1/2.04;

	fluxes.gpp    = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.gpp;});
	fluxes.npp    = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.npp;});
	fluxes.trans  = m1      * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.trans;});
	fluxes.tleaf  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.tleaf;});
	fluxes.troot  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.troot;});
	fluxes.rleaf  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rleaf;});
	fluxes.rroot  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rroot;});
	fluxes.rstem  = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->res.rstem;});
	fluxes.mort   = m1 * m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->get_biomass()*p->rates.dmort_dt;});
	fluxes.gs = (fluxes.trans*55.55/86400)/1.6/(static_cast<PSPM_Environment*>(S.env)->clim_inst.vpd/1.0325e5);
	//     ^ convert kg/m2/day --> mol/m2/s

	// Note: for these vector calcs, multipliers are multiplied inside lambdas, where necesary
	species.n_ind_vec       = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant * p){return 1;});
	species.biomass_vec     = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant * p){return p->get_biomass();});
	species.basal_area_vec  = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant * p){
		                                                                        double D = p->geometry.diameter_at_height(1.3, p->traits);
	                                                                            return M_PI*D*D/4;
	                                                                         });
	species.canopy_area_vec = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant * p){return p->geometry.crown_area;});
	species.height_vec      = integrate_prop_above_per_species(t, 0.1, S, [](PSPM_Plant * p){return p->geometry.height;});
	for (int k=0; k<S.n_species(); ++k) species.height_vec[k] /= species.n_ind_vec[k];
	species.mortality_vec   = integrate_prop_per_species(t, S, [m1, m2](PSPM_Plant * p){return m1*m2*p->get_biomass()*p->rates.dmort_dt;});

	structure.n_ind         = std::accumulate(species.n_ind_vec.begin(),       species.n_ind_vec.end(),       0.0);
	structure.biomass       = std::accumulate(species.biomass_vec.begin(),     species.biomass_vec.end(),     0.0);
	structure.basal_area    = std::accumulate(species.basal_area_vec.begin(),  species.basal_area_vec.end(),  0.0);
	structure.canopy_area   = std::accumulate(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), 0.0);

	structure.lai           =      integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.crown_area*p->geometry.lai;});
	structure.leaf_mass     = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.leaf_mass(p->traits);});
	structure.stem_mass     = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.stem_mass(p->traits);});
	structure.croot_mass    = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.coarse_root_mass(p->traits);});
	structure.froot_mass    = m2 * integrate_prop(t, S, [](PSPM_Plant* p){return p->geometry.root_mass(p->traits);});

	// upper canopy acclimated traits
	double h_uc = static_cast<PSPM_Environment*>(S.env)->z_star[0];
	structure.canopy_area_uc = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant * p){return (p->geometry.height > h_uc)? p->geometry.crown_area : 0;});

	acc_traits.vcmax = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant * p){return (p->geometry.height > h_uc)? (p->res.vcmax25_avg * p->geometry.crown_area) : 0; });
	acc_traits.vcmax /= structure.canopy_area_uc;

	acc_traits.dpsi  = integrate_prop_above(t, 0.1, S, [h_uc](PSPM_Plant * p){return (p->geometry.height > h_uc)? (p->res.dpsi_avg * p->geometry.crown_area) : 0; });
	acc_traits.dpsi /= structure.canopy_area_uc;

	misc.cc_est = (fluxes.tleaf + fluxes.troot + fluxes.rleaf + fluxes.rroot + fluxes.rstem)/fluxes.tleaf;

	// LAI vertical profile
	structure.lai_vert.clear();
	structure.lai_vert.resize(25, 0);
	for (int iz=0; iz<25; ++iz){
		structure.lai_vert[iz] = integrate_prop(t, S, [iz](PSPM_Plant * p){return p->geometry.crown_area_above(iz,p->traits)*p->geometry.lai;});
	}
}


CommunityProperties operator + (CommunityProperties lhs, CommunityProperties &rhs){
	lhs += rhs;
	return lhs;
}


void CommunityProperties::openStreams(std::string dir){

	if (b_output_cohort_props){
		cohort_props_out.open(dir + "/cohort_props.txt");
		cohort_props_out << "t\tspeciesID\tcohortID\t";
		for (auto vname : varnames) cohort_props_out << vname << "\t";
		cohort_props_out << std::endl;
		cohort_props_out << std::setprecision(12);
	}

	size_dists_out.open(dir + "/size_distributions.txt");
	fzst.open(std::string(dir + "/z_star.txt").c_str());
	fco.open(std::string(dir + "/canopy_openness.txt").c_str());
	flai.open(std::string(dir + "/lai_profile.txt").c_str());
	foutd.open(std::string(dir + "/" + emgProps_file).c_str());
	fouty.open(std::string(dir + "/" + cwmAvg_file).c_str());
	fouty_spp.open(std::string(dir + "/" + cwmperSpecies_file).c_str());
	ftraits.open(std::string(dir + "/" + traits_file).c_str());
	// fclim.open(std::string(dir + "/climate_co2.txt").c_str());

	foutd << "YEAR\tGPP\tNPP\tRAU\tMORT\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\tVCMAX\tDPSI\tCCEST\tCO2\n";
	fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
	fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\tSEEDS\n";
	ftraits << "YEAR\tSPP\tRES\tLMA\tWD\tHMAT\tP50X\tZETA\tr0_last\tr0_avg\tr0_exp\tr0_cesaro\n";
	// fclim << "t\ttc\tppfd_max\tppfd\tvpd\tco2\telv\tswp\n";

}

void CommunityProperties::closeStreams(){
	if (b_output_cohort_props) cohort_props_out.close();
	size_dists_out.close();
	
	fzst.close();
	fco.close();
	flai.close();
	foutd.close();
	fouty.close();
	fouty_spp.close();
	ftraits.close();
	// fclim.close();
}

void CommunityProperties::writeOut(double t, Patch &P){
	Solver *S = &P.S;

	// consistently output date in decimal years across all files
	double date = flare::julian_to_yearsCE(P.ts.to_julian(t));

	for (int s=0; s < S->species_vec.size(); ++s){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[s]);

		std::vector<double> breaks = my_log_seq(0.01, 10, 100);
		std::vector<double> dist = S->getDensitySpecies1D(s, 0, breaks);
		//cout << "here: " << breaks.size() << " " << dist.size() << endl;

		if (spp->isResident){
			size_dists_out << date << "\t" << spp->species_name << "\t";
			for (int i=0; i<100; ++i){
				size_dists_out << dist[i] << "\t";
			}
			size_dists_out << "\n";
		}
		
		if (b_output_cohort_props){
			if (spp->isResident){
				for (int j=0; j<spp->xsize()-1; ++j){
					auto& C = spp->getCohort(j);
					cohort_props_out 
						<< date << "\t" 
						<< spp->species_name << "\t"  // use name instead of index s becuase it is unique and order-insensitive
						<< j << "\t"
						<< C.geometry.diameter << "\t"
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
	}

	foutd 
		<< date << "\t"
		<< fluxes.gpp << "\t"
		<< fluxes.npp << "\t"
		<< (fluxes.rleaf + fluxes.rroot + fluxes.rstem) << "\t"  // kgC/m2/d
		<< fluxes.mort << "\t"
		<< structure.leaf_mass << "\t"     
		<< structure.stem_mass << "\t"
		<< structure.croot_mass << "\t"
		<< structure.froot_mass << "\t"
		<< (structure.croot_mass+structure.froot_mass) << "\t" // kgC/m2
		<< fluxes.gs << "\t"
		<< fluxes.trans << "\t"   
		<< structure.lai << "\t"
		<< acc_traits.vcmax << "\t"
		<< acc_traits.dpsi << "\t"
		<< misc.cc_est << "\t"
		<< static_cast<PSPM_Environment*>(S->env)->clim_inst.co2 
		<< std::endl;
	
	fouty 
		<< date << "\t"
		<< -9999  << "\t"
		<< structure.n_ind << "\t"
		<< -9999  << "\t"
		<< -9999  << "\t" // height
		<< -9999 /*cwm.hmat*/  << "\t"
		<< structure.canopy_area  << "\t"   // m2/m2
		<< structure.basal_area  << "\t"            // m2/m2
		<< structure.biomass  << "\t"       // kg/m2
		<< -9999 /*cwm.wd*/  << "\t"
		<< -9999  << "\t"
		<< -9999 /*1/cwm.lma*/  << "\t"
		<< -9999 /*cwm.p50*/  
		<< std::endl;
	
	for (int k=0; k<S->species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[k]);
		fouty_spp 
			<< date << "\t"
			<< spp->species_name  << "\t" // use name instead of index k becuase it is unique and order-insensitive
			<< species.n_ind_vec[k] << "\t"
			<< -9999  << "\t"
			<< species.height_vec[k]  << "\t"
			<< -9999 /*species.hmat_vec[k]*/  << "\t"
			<< species.canopy_area_vec[k]  << "\t"   // m2/m2
			<< species.basal_area_vec[k]  << "\t"            // m2/m2
			<< species.biomass_vec[k]  << "\t"       // kg/m2
			<< -9999 /*species.wd_vec[k]*/  << "\t"
			<< -9999  << "\t"
			<< -9999 /*1/cwm.lma_vec[k]*/  << "\t"
			<< -9999 /*cwm.p50_vec[k]*/  << "\t"
			<< spp->seeds_hist.get()  
			<< std::endl;
	}

	for (int k=0; k<S->species_vec.size(); ++k){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[k]);
		ftraits 
			<< date << "\t"
			<< spp->species_name  << "\t" // use name instead of index k becuase it is unique and order-insensitive
			<< spp->isResident << "\t"
			<< spp->getCohort(-1).traits.lma << "\t"
			<< spp->getCohort(-1).traits.wood_density << "\t"
			<< spp->getCohort(-1).traits.hmat << "\t"
		    << spp->getCohort(-1).traits.p50_xylem << "\t"
			<< spp->getCohort(-1).traits.zeta << "\t"
			<< spp->r0_hist.get_last() << "\t"
			<< spp->r0_hist.get() << "\t"
			<< 0 << "\t"
			<< 0 //spp->r0_hist.get_cesaro();
			<< std::endl; 
	}

	ftraits.flush();
	fouty_spp.flush();

	flai << date << "\t";
	for (int i=0; i<structure.lai_vert.size(); ++i) flai << structure.lai_vert[i] << "\t";
	flai << std::endl;	

	fzst << date << "\t";
	for (auto z : static_cast<PSPM_Environment*>(S->env)->z_star) fzst << z << "\t";
	fzst << std::endl;
	
	fco << date << "\t";
	for (auto z : static_cast<PSPM_Environment*>(S->env)->canopy_openness) fco << z << "\t";
	fco << std::endl;

}

// FIXME: These operators are not tested
CommunityProperties& CommunityProperties::operator /= (double s){

	fluxes.gpp            /= s;
	fluxes.npp            /= s;
	fluxes.trans          /= s;
	fluxes.gs             /= s;
	fluxes.tleaf          /= s;
	fluxes.troot          /= s;
	fluxes.rleaf          /= s;
	fluxes.rroot          /= s;
	fluxes.rstem          /= s;

	structure.leaf_mass   /= s;
	structure.stem_mass   /= s;
	structure.croot_mass  /= s;
	structure.froot_mass  /= s;
	structure.biomass     /= s;
	structure.basal_area  /= s;
	structure.canopy_area /= s;
	structure.n_ind       /= s;
	structure.height      /= s;
	structure.lai         /= s;
	transform(structure.lai_vert.begin(), structure.lai_vert.end(), structure.lai_vert.begin(), [s](const double &c){ return c/s; });

	misc.cc_est           /= s;

	acc_traits.vcmax      /= s;
	acc_traits.dpsi       /= s;

	transform(species.n_ind_vec.begin(),       species.n_ind_vec.end(),       species.n_ind_vec.begin(),       [s](const double &c){ return c/s; });
	transform(species.biomass_vec.begin(),     species.biomass_vec.end(),     species.biomass_vec.begin(),     [s](const double &c){ return c/s; });
	transform(species.basal_area_vec.begin(),  species.basal_area_vec.end(),  species.basal_area_vec.begin(),  [s](const double &c){ return c/s; });
	transform(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), species.canopy_area_vec.begin(), [s](const double &c){ return c/s; });
	transform(species.height_vec.begin(),      species.height_vec.end(),      species.height_vec.begin(),      [s](const double &c){ return c/s; });
	
	return *this;
}

CommunityProperties& CommunityProperties::operator += (const CommunityProperties &rhs){

	fluxes.gpp            += rhs.fluxes.gpp;
	fluxes.npp            += rhs.fluxes.npp;
	fluxes.trans          += rhs.fluxes.trans;
	fluxes.gs             += rhs.fluxes.gs;
	fluxes.tleaf          += rhs.fluxes.tleaf;
	fluxes.troot          += rhs.fluxes.troot;
	fluxes.rleaf          += rhs.fluxes.rleaf;
	fluxes.rroot          += rhs.fluxes.rroot;
	fluxes.rstem          += rhs.fluxes.rstem;

	structure.leaf_mass   += rhs.structure.leaf_mass;
	structure.stem_mass   += rhs.structure.stem_mass;
	structure.croot_mass  += rhs.structure.croot_mass;
	structure.froot_mass  += rhs.structure.froot_mass;
	structure.biomass     += rhs.structure.biomass;
	structure.basal_area  += rhs.structure.basal_area;
	structure.canopy_area += rhs.structure.canopy_area;
	structure.n_ind       += rhs.structure.n_ind;
	structure.height      += rhs.structure.height;
	structure.lai         += rhs.structure.lai;
	transform(structure.lai_vert.begin(), structure.lai_vert.end(), rhs.structure.lai_vert.begin(), structure.lai_vert.begin(), std::plus<double>());

	misc.cc_est           += rhs.misc.cc_est;

	acc_traits.vcmax      += rhs.acc_traits.vcmax;
	acc_traits.dpsi       += rhs.acc_traits.dpsi;

	transform(species.n_ind_vec.begin(),       species.n_ind_vec.end(),       rhs.species.n_ind_vec.begin(),       species.n_ind_vec.begin(),       std::plus<double>());
	transform(species.biomass_vec.begin(),     species.biomass_vec.end(),     rhs.species.biomass_vec.begin(),     species.biomass_vec.begin(),     std::plus<double>());
	transform(species.basal_area_vec.begin(),  species.basal_area_vec.end(),  rhs.species.basal_area_vec.begin(),  species.basal_area_vec.begin(),  std::plus<double>());
	transform(species.canopy_area_vec.begin(), species.canopy_area_vec.end(), rhs.species.canopy_area_vec.begin(), species.canopy_area_vec.begin(), std::plus<double>());
	transform(species.height_vec.begin(),      species.height_vec.end(),      rhs.species.height_vec.begin(),      species.height_vec.begin(),      std::plus<double>());
	
	return *this;
}



} // namespace pfate