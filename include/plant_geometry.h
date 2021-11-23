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
	private:
	struct{
		// geometry traits
		double m, n;		// crown shape paramaters
		double a;			// height-diameter allometry
		double c;			// crown area allometry
		double fg;			// upper canopy gap fraction
		
		// Precomputed Geometric parameters
		double eta_c;
		double pic_4a;
		double zm_H;
		double qm;
	} geom;

	public:
	// current state
	double height;	// height
	double diameter;	// basal diameter
	double crown_area;	// crown area
	double leaf_area;	// leaf area
	double sapwood_fraction = 1;	// sapwood fraction
	double r0;

	//double fl = 1; // current realized fraction of the maximum possible LAI
	//double hvlc = 1;

	double sap_frac_ode = 1;
	double heart_mass_ode = 0;
	double k_sap;

	public:

	void initGeometry(double a, double c, double m, double n, double fg){
		geom.m = m; geom.n = n; 
		geom.a = a; geom.c = c;
		geom.fg = fg;

		geom.pic_4a = M_PI*geom.c/(4*geom.a);

		geom.zm_H = pow((n-1)/(m*n-1), 1/n);
		geom.qm = m*n * pow((n-1)/(m*n-1), 1-1/n) * pow((m-1)*n/(m*n-1), m-1);

		geom.eta_c = geom.zm_H - m*m*n/(geom.qm*geom.qm) * beta(2-1/n, 2*m-1) * (incbeta(2-1/n, 2*m-1, (n-1)/(m*n-1)) - (1-geom.fg)); 
		
		std::cout << "m = " << m << ", n = " << n << ", zm/H = " << geom.zm_H << ", qm = " << geom.qm << ", eta_c = " << geom.eta_c << "\n";
	}


	//void set_size(double _x, PlantTraits &traits){
		//height = _x;
		//double lai = geom.lai_max * traits.fl;
		//double hv_min = 1/(par.lai_max * par.c);
		//double hv = 1/(lai*par.c);
		//diameter = -log(1-height/traits.hmat) * traits.hmat/par.a;
		//crown_area = par.pic_4a * height * diameter;
		//leaf_area = crown_area*lai;
		//sapwood_fraction = (hvlc) * height/diameter/par.a;	// FIXME: check LAI variation
	//}

	//double dsize_dmass(PlantTraits &traits) const {
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

	void set_size(double _x, PlantTraits &traits){
		diameter = _x;
		height = traits.hmat * (1 - exp(-geom.a*diameter/traits.hmat));
		crown_area = geom.pic_4a * height * diameter;
		leaf_area = crown_area * traits.lai;
		sapwood_fraction = height / (diameter * geom.a);	
		r0 = sqrt(crown_area/M_PI)/geom.qm; 
	}


	double dsize_dmass(PlantTraits &traits) const {
		double dh_dd = geom.a * exp(-geom.a*diameter/traits.hmat);
		double dmleaf_dd = traits.lma * traits.lai * geom.pic_4a * (height + diameter*dh_dd);	// FIXME: Carefully check LAI variation
		double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
		double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
		double dmroot_dd = (traits.zeta/traits.lma) * dmleaf_dd;

		double dmass_dd = dmleaf_dd + dmtrunk_dd + dmbranches_dd + dmroot_dd;
		return 1/dmass_dd;
	}


	double leaf_mass(PlantTraits &traits){
		return leaf_area*traits.lma;	
	}

	double root_mass(PlantTraits &traits){
		return leaf_area*traits.zeta;	
	}
	
	//double sapwood_mass(PlantTraits &traits){
		//return traits.wood_density*(hvlc/geom.c)*crown_area*geom.eta_l*height;
	//}
		
	double sapwood_mass(PlantTraits &traits){
		return stem_mass(traits)*sapwood_fraction;
	}

	double stem_mass(PlantTraits &traits){
		double trunk_mass = traits.wood_density*(M_PI*diameter*diameter/4)*height*geom.eta_c;
		double branch_mass = traits.wood_density * (M_PI*diameter*diameter/12)*height * sqrt((geom.c/geom.a)*(diameter/height));	
		return trunk_mass + branch_mass;
	}

	double heartwood_mass(PlantTraits &traits){
		return stem_mass(traits)*(1-sapwood_fraction);
	}

	double total_mass(PlantTraits &traits){
		return stem_mass(traits) + leaf_mass(traits) + root_mass(traits);
	}


	double q(double z, PlantParameters &par){
		if (z > height || z < 0) return 0;
		else{
			double m = geom.m, n = geom.n;
			double zHn_1 = pow(z/height, n-1);
			double zHn   = zHn_1 * z/height;
			return m*n * pow(1 - zHn, m-1) * zHn_1;
		}
	}

	double zm(PlantParameters &par){
		return geom.zm_H * height;
	} 

	//double Q(double z, double H, double n, double m){
		//if (z > H) return 0;
		//else if (z <= 0) return 1;
		//else{
			//return pow(1-pow(z/H, n), m); 
		//}
	//}

	
	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - simulates growth over dt with constant assimilation rate A
	// ** 
	void grow_for_dt(double t, double dt, double &prod, double A, PlantTraits &traits){

		auto derivs = [A, &traits, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
			set_size(S[1], traits);

			double dh_dd = geom.a*exp(-geom.a*diameter/traits.hmat);

			dSdt[0] = A*leaf_area;	// biomass production rate
			dSdt[1] = dsize_dmass(traits) * A*leaf_area; 
			dSdt[2] = 1/(geom.a*diameter*diameter)*(diameter*dh_dd - height)*dSdt[1];
			
			double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
			double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
			
			double dsap_trunk_dd = traits.wood_density * M_PI / (4*geom.a) * geom.eta_c * (2*diameter*dh_dd + height) * height;
			double dsap_branch_dd = traits.wood_density * M_PI / (8*geom.a) * sqrt(geom.c/geom.a) * (diameter*dh_dd + height) * sqrt(diameter*height);

			dSdt[3] = (dmtrunk_dd + dmbranches_dd - dsap_trunk_dd - dsap_branch_dd) * dSdt[1];
				
			k_sap = dSdt[3]/sapwood_mass(traits)*dSdt[1];
		};

		std::vector<double> S = {prod, get_size(), sap_frac_ode, heart_mass_ode};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		heart_mass_ode = S[3];
		sap_frac_ode = S[2];
		set_size(S[1], traits);
		prod = S[0];
	}


};


} // namespace plant

#endif


