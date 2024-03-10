#ifndef PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_
#define PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_

#include "pspm_interface.h"
#include "adaptive_species.h"
#include "utils/sequence.h"

#ifndef M_PI
#define M_PI 3.14159265358
#endif

namespace pfate {

class Patch;

class CommunityProperties{
	public:
	struct{
		double gpp=0;
		double npp=0;
		double trans=0;
		double gs=0;
		double tleaf=0;
		double troot=0;
		double rleaf=0;
		double rroot=0;
		double rstem=0;
	} fluxes;

	struct{
		double leaf_mass=0;
		double stem_mass=0;
		double croot_mass=0;
		double froot_mass=0;
		double biomass=0;
		double basal_area=0;
		double canopy_area=0;
		double canopy_area_uc=0;
		double n_ind=0;
		double height=0;
		double lai=0;
		std::vector <double> lai_vert;
	} structure;

	struct{
		std::vector<double> n_ind_vec;
		std::vector<double> biomass_vec;
		std::vector<double> basal_area_vec;
		std::vector<double> canopy_area_vec;
		std::vector<double> height_vec;
	} species;

	struct{
		double cc_est=0;
	} misc;

	struct{
		double vcmax=0;
		double dpsi=0;
	} acc_traits;

	// double lma=0;
	// double p50=0;
	// double hmat=0;
	// double wd=0;
	
	// std::vector<double> lma_vec;
	// std::vector<double> p50_vec;
	// std::vector<double> hmat_vec;
	// std::vector<double> wd_vec;
	
	bool b_output_cohort_props = false;
	std::string emgProps_file;
	std::string cwmAvg_file;
	std::string cwmperSpecies_file;
	std::string traits_file;

	public:
	void openStreams(std::string dir);
	void closeStreams();
	void writeOut(double t, Patch &P);

	void resize(int n);
	bool isResident(Species_Base * spp);
	void update(double t, Patch &P);
	
	CommunityProperties & operator /= (double s);
	CommunityProperties & operator += (const CommunityProperties &s);

	private:
	std::vector<std::string> varnames = {"diameter", "height", "lai", "mort", "fec", "rgr", "gpp"};
	std::ofstream cohort_props_out;
	std::ofstream size_dists_out;
	std::ofstream fzst;
	std::ofstream fco;
	std::ofstream flai;
	std::ofstream foutd;
	std::ofstream fouty;
	std::ofstream fouty_spp;
	std::ofstream ftraits;
	std::ofstream fclim;


	template<class Func>
	double integrate_prop(double t, Solver &S, Func f){
		double x = 0;
		for (int k=0; k<S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if (is_resident){
				x += S.state_integral([&S,k,f](int i, double t){
						auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
						PSPM_Plant * pp = &p;
						return f(pp);
					}, t, k);
			}
		}
		return x;
	}

	template<class Func>
	std::vector<double> integrate_prop_per_species(double t, Solver &S, Func f){
		std::vector<double> x; 
		x.reserve(S.n_species());
		for (int k=0; k<S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if (is_resident){
				x.push_back(
					S.state_integral([&S,k,f](int i, double t){
						auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
						PSPM_Plant * pp = &p;
						return f(pp);
					}, t, k)
				);
			}
		}
		return x;
	}

	template<class Func>
	std::vector<double> integrate_prop_above_per_species(double t, double xlow, Solver &S, Func f){
		std::vector<double> x; 
		x.reserve(S.n_species());
		for (int k=0; k<S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if (is_resident){
				x.push_back(
					S.integrate_wudx_above([&S,k,f](int i, double t){
						auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
						PSPM_Plant * pp = &p;
						return f(pp);
					}, t, {xlow}, k)
				);
			}
		}
		return x;
	}

	template<class Func>
	double integrate_prop_above(double t, double xlow, Solver &S, Func f){
		double x = 0;
		for (int k=0; k<S.n_species(); ++k){
			bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
			if (is_resident){
				x += S.integrate_wudx_above([&S,k,f](int i, double t){
						auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
						PSPM_Plant * pp = &p;
						return f(pp);
					}, t, {xlow}, k);
			}
		}
		return x;
	}


};

CommunityProperties operator + (CommunityProperties lhs, CommunityProperties &rhs);

} // namespace pfate

#endif
