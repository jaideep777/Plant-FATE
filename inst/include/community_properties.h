#ifndef PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_
#define PLANT_FATE_PFATE_COMMUNITY_PROPERTIES_H_

#include "pspm_interface.h"
#include "adaptive_species.h"
#include "utils/sequence.h"

#ifndef M_PI
#define M_PI 3.14159265358
#endif

namespace pfate {

template<class Func>
double integrate_prop(double t, Solver &S, const Func &f){
	double x = 0;
	for (int k=0; k<S.n_species(); ++k){
		bool is_resident = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S.species_vec[k])->isResident;
		if (is_resident){
			x += S.state_integral([&S,k,f](int i, double t){
					auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
					const PSPM_Plant * pp = &p;
					return f(pp);
				 }, t, k);
		}
	}
	return x;
}



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
	double dpsi=0;

	std::vector<double> n_ind_vec;
	std::vector<double> biomass_vec;
	std::vector<double> ba_vec;
	std::vector<double> canopy_area_vec;
	std::vector<double> height_vec;
	std::vector<double> vcmax_vec;
	std::vector<double> dpsi_vec;
	
	std::vector<double> lma_vec;
	std::vector<double> p50_vec;
	std::vector<double> hmat_vec;
	std::vector<double> wd_vec;
	
	void resize(int n);
	bool isResident(Species_Base * spp);
	void update(double t, Solver &S);
	
	SpeciesProps & operator /= (double s);
	SpeciesProps & operator += (const SpeciesProps &s);

};

SpeciesProps operator + (SpeciesProps lhs, SpeciesProps &rhs);


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
	double cc_est=0;

	std::vector <double> lai_vert;

	
	EmergentProps & operator /= (double s);
	
	EmergentProps & operator += (const EmergentProps &s);

	bool isResident(Species_Base * spp);
	
	void update(double t, Solver &S);

};

EmergentProps operator + (EmergentProps lhs, EmergentProps &rhs);


class SolverIO{
	public:
	int nspecies;
	Solver * S;
	std::vector<std::string> varnames = {"diameter", "height", "lai", "mort", "fec", "rgr", "gpp"};

	bool b_output_cohort_props = false;

	std::string emgProps_file;
	std::string cwmAvg_file;
	std::string cwmperSpecies_file;
	std::string traits_file;

	// std::vector <std::vector<std::ofstream>> streams;
	std::ofstream cohort_props_out;
	std::ofstream size_dists_out;

	std::ofstream fzst;
	std::ofstream fco;
	// std::ofstream fseed;   // not using these 2 because they are not insensitive to order of species
	// std::ofstream fabase;
	std::ofstream flai;
	std::ofstream foutd;
	std::ofstream fouty;
	std::ofstream fouty_spp;
	std::ofstream ftraits;

	std::ofstream fclim;

	void openStreams(std::string dir);

	void closeStreams();

	void writeState(double t, SpeciesProps& cwm, EmergentProps& props);
};

} // namespace pfate

#endif
