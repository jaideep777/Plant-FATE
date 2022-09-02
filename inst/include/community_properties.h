#ifndef PLANT_FATE_COMMUNITY_PROPERTIES_H_
#define PLANT_FATE_COMMUNITY_PROPERTIES_H_

#include <solver.h>

#include <solver.h>

// FIXME: move definitions to cpp
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
	

	void update(double t, Solver &S){
		n_ind_vec.clear();
		n_ind_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			n_ind_vec[k] = S.integrate_x([&S,k](int i, double t){
										      return 1;
										}, t, k);
		n_ind = std::accumulate(n_ind_vec.begin(), n_ind_vec.end(), 0.0);

		biomass_vec.clear();
		biomass_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			biomass_vec[k] = S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.get_biomass();
										}, t, k);
		biomass = std::accumulate(biomass_vec.begin(), biomass_vec.end(), 0.0);

		ba_vec.clear();
		ba_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			ba_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
											  double D = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry.diameter;
											  return M_PI*D*D/4;
										}, t, 0.1, k);
		ba = std::accumulate(ba_vec.begin(), ba_vec.end(), 0.0);

		canopy_area_vec.clear();
		canopy_area_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			canopy_area_vec[k] = S.integrate_x([&S,k](int i, double t){
											  auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry;
											  return p.crown_area;
										}, t, k);
		canopy_area = std::accumulate(canopy_area_vec.begin(), canopy_area_vec.end(), 0.0);
		

		height_vec.clear();
		height_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			height_vec[k] = S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.geometry.height;
										}, t, k);
										
		for (int k=0; k<S.n_species(); ++k) height_vec[k] /= n_ind_vec[k];


		hmat = 0;
		for (int k=0; k<S.n_species(); ++k)
			hmat += S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.traits.hmat;
										}, t, k);
		hmat /= n_ind;
		hmat_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k) hmat_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.hmat;


		lma = 0;
		for (int k=0; k<S.n_species(); ++k)
			lma += S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.traits.lma;
										}, t, k);
		lma /= n_ind;
		lma_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k) lma_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.lma;

		wd = 0;
		for (int k=0; k<S.n_species(); ++k)
			wd += S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.traits.wood_density;
										}, t, k);
		wd /= n_ind;
		wd_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k) wd_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.wood_density;

		p50 = 0;
		for (int k=0; k<S.n_species(); ++k)
			p50 += S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.traits.p50_xylem;
										}, t, k);
		p50 /= n_ind;
		p50_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k) p50_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.p50_xylem;

		gs = 0;
		for (int k=0; k<S.n_species(); ++k)
			gs += S.integrate_x([&S,k](int i, double t){
										      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
										      return p.res.gs_avg * p.geometry.crown_area;
										}, t, k);
		gs /= canopy_area;

		vcmax_vec.clear();
		vcmax_vec.resize(S.n_species());
		for (int k=0; k<S.n_species(); ++k)
			vcmax_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
											  auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
											  return p.res.vcmax_avg * p.geometry.crown_area;
										}, t, 0.1, k);
		vcmax = std::accumulate(vcmax_vec.begin(), vcmax_vec.end(), 0.0);
		vcmax /= canopy_area;
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
	
	vector <double> lai_vert;

	
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
		
		transform(lai_vert.begin(), lai_vert.end(), lai_vert.begin(), [s](const double &c){ return c/s; });
		
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

		transform(lai_vert.begin(), lai_vert.end(), s.lai_vert.begin(), lai_vert.begin(), std::plus<double>());

		return *this;
	}
	
	template<class Func>
	double integrate_prop(double t, Solver &S, const Func &f){
		double x = 0;
		for (int k=0; k<S.n_species(); ++k)
			x += S.integrate_x([&S,k,f](int i, double t){
				auto p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
				const PSPM_Plant * pp = &p;
				return f(pp);
			}, t, k);
		return x;
	}
	
	void update(double t, Solver &S){
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

		// LAI vertical profile
		lai_vert.resize(25);
		for (int iz=0; iz<25; ++iz)
			for (int k=0; k<S.n_species(); ++k)
				lai_vert[iz] += S.integrate_x([&S,k,iz](int i, double t){
									  auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
									  return p.geometry.crown_area_above(iz,p.traits)*p.geometry.lai;
								}, t, k);

	}

};

EmergentProps operator + (EmergentProps lhs, EmergentProps &rhs){
	lhs += rhs;
	return lhs;
}


#endif