#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#include <cmath>

#include "plant_params.h"

namespace plant{

template <class functor, class container>
void Euler(double x, double h, container& y, functor& derivs){
	container fk(y.size());
	derivs(x, y, fk);
	for (int i=0; i<y.size(); i++) y[i] += h*fk[i]; 
}

template <class functor, class container>
void RK4(double x, double h, container& y, functor& derivs){
	static container k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());// temporary arrays
	static container yt(y.size());

	double h2=h*0.5;
	double xh = x + h2;
	derivs(x, y, k1);     // First step : evaluating k1
	for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k1[i];// Preparing second step by  ty <- y + k1/2
	derivs(xh, yt, k2);                                    // Second step : evaluating k2
	for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k2[i];// Preparing third step by   yt <- y + k2/2
	derivs(xh, yt, k3);                                    // Third step : evaluating k3
	for (int i=0; i<y.size(); i++) yt[i] = y[i] +  h*k3[i];// Preparing fourth step  yt <- y + k3
	derivs(x+h, yt, k4);                                   // Final step : evaluating k4
	for (int i=0; i<y.size(); i++) y[i] += h/6.0*(k1[i]+2.0*(k2[i]+k3[i])+k4[i]);
}


class PlantGeometry{
	public:
	double height;	// height
	double diameter;	// basal diameter
	double crown_area;	// crown area
	double leaf_area;	// leaf area
	double sapwood_fraction = 1;	// sapwood fraction

	//double fl = 1; // current realized fraction of the maximum possible LAI
	double hvlc = 1;

	//void set_size(double _x, PlantParameters &par, PlantTraits &traits){
		//height = _x;
		//double lai = par.lai_max * traits.fl;
		//double hv_min = 1/(par.lai_max * par.c);
		//double hv = 1/(lai*par.c);
		//diameter = -log(1-height/traits.hmat) * traits.hmat/par.a;
		//crown_area = par.pic_4a * height * diameter;
		//leaf_area = crown_area*lai;
		//sapwood_fraction = (hvlc) * height/diameter/par.a;	// FIXME: check LAI variation
	//}

	//double dsize_dmass(PlantParameters &par, PlantTraits &traits) const {
		//double lai = par.lai_max * traits.fl;
		//double dd_dh = 1/(par.a*(1-height/traits.hmat));
		//double dmleaf_dh = traits.lma*lai*par.pic_4a * (height*dd_dh + diameter);	// FIXME: Carefully check LAI variation
		//double dmstem_dh = (par.eta_l*M_PI*traits.wood_density/4) * (2*height*dd_dh + diameter)*diameter;
		//double dmroot_dh = (traits.zeta/traits.lma) * dmleaf_dh;

		//double dmass_dh = dmleaf_dh + dmstem_dh + dmroot_dh;
		//return 1/dmass_dh;
	//}


	double get_size() const {
		return diameter;
	}

	void set_size(double _x, PlantParameters &par, PlantTraits &traits){
		diameter = _x;
		double lai = par.lai_max * traits.fl;
		height = traits.hmat * (1 - exp(-par.a*diameter/traits.hmat));
		crown_area = par.pic_4a * height * diameter;
		leaf_area = crown_area*lai;
		sapwood_fraction = (hvlc) * height/diameter/par.a;	// FIXME: check LAI variation
	}


	double dsize_dmass(PlantParameters &par, PlantTraits &traits) const {
		double lai = par.lai_max * traits.fl;
		double dh_dd = par.a*exp(-par.a*diameter/traits.hmat);
		double dmleaf_dd = traits.lma*lai*par.pic_4a * (height + diameter*dh_dd);	// FIXME: Carefully check LAI variation
		double dmstem_dd = (par.eta_l*M_PI*traits.wood_density/4) * (2*height + diameter*dh_dd)*diameter;
		double dmroot_dd = (traits.zeta/traits.lma) * dmleaf_dd;

		double dmass_dd = dmleaf_dd + dmstem_dd + dmroot_dd;
		return 1/dmass_dd;
	}


	double leaf_mass(PlantParameters &par, PlantTraits &traits){
		return leaf_area*traits.lma;	
	}

	double root_mass(PlantParameters &par, PlantTraits &traits){
		return leaf_area*traits.zeta;	
	}
	
	double sapwood_mass(PlantParameters &par, PlantTraits &traits){
		return traits.wood_density*(hvlc/par.c)*crown_area*par.eta_l*height;
	}
		
	double sapwood_mass1(PlantParameters &par, PlantTraits &traits){
		return stem_mass(par, traits)*sapwood_fraction;
	}

	double stem_mass(PlantParameters &par, PlantTraits &traits){
		return traits.wood_density*(M_PI*diameter*diameter/4)*height*par.eta_l;	
	}

	double total_mass(PlantParameters &par, PlantTraits &traits){
		return stem_mass(par, traits) + leaf_mass(par, traits) + root_mass(par, traits);
	}


	double q(double z, double H, double n, double m){
		if (z > H || z < 0) return 0;
		else if (z == 0){
			if (n == 1) return m/H;
			else return 0;
		}
		else{
			double zH_n = pow(z/H, n);
			return m*n*pow(1-zH_n, m-1)*zH_n/z;
		}
	}

	double Q(double z, double H, double n, double m){
		if (z > H) return 0;
		else if (z <= 0) return 1;
		else{
			return pow(1-pow(z/H, n), m); 
		}
	}

	
	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - simulates growth over dt with constant assimilation rate A
	// ** 
	void grow_for_dt(double t, double dt, double &prod, double A, PlantParameters &par, PlantTraits &traits){

		auto derivs = [A, &par, &traits, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			set_size(S[1], par, traits);

			dSdt[0] = A*leaf_area;	// biomass production rate
			dSdt[1] = dsize_dmass(par, traits) * A*leaf_area; 
		};

		std::vector<double> S = {prod, get_size()};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		set_size(S[1], par, traits);
		prod = S[0];
	}


};


} // namespace plant

#endif


