#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#include <cmath>

#include "lambertw.h"
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
	double lai = 1;     // leaf area index 
	double diameter;	// basal diameter
	double litter_pool = 0;

	double height;	// height
	double crown_area;	// crown area
	//double leaf_area;	// leaf area
	double sapwood_fraction = 1;	// fraction of stem that is sapwood
	double functional_xylem_fraction;	// fraction of funcitonal xylem in sapwood

	double sap_frac_ode = 1;
	double sapwood_mass_ode = 0;
	double heart_mass_ode = 0;
	double k_sap;

	public:

	void initGeometry(double diameter_0, PlantParameters &par, PlantTraits &traits){
		geom.m = par.m; geom.n = par.n; 
		geom.a = par.a; geom.c = par.c;
		geom.fg = par.fg;

		geom.pic_4a = M_PI*geom.c/(4*geom.a);

		double m = geom.m, n = geom.n;
		geom.zm_H = pow((n-1)/(m*n-1), 1/n);
		geom.qm = m*n * pow((n-1)/(m*n-1), 1-1/n) * pow((m-1)*n/(m*n-1), m-1);

		geom.eta_c = geom.zm_H - m*m*n/(geom.qm*geom.qm) * beta(2-1/n, 2*m-1) * (incbeta(2-1/n, 2*m-1, (n-1)/(m*n-1)) - (1-geom.fg)); 
		
		std::cout << "m = " << m << ", n = " << n << ", zm/H = " << geom.zm_H << ", qm = " << geom.qm << ", eta_c = " << geom.eta_c << "\n";
		
		set_size(diameter_0, traits);
	}

	double q(double z){
		if (z > height || z < 0) return 0;
		else{
			double m = geom.m, n = geom.n;
			double zHn_1 = pow(z/height, n-1);
			double zHn   = zHn_1 * z/height;
			return m*n * pow(1 - zHn, m-1) * zHn_1;
		}
	}

	double zm(){
		return geom.zm_H * height;
	} 

	
	double get_size() const {
		return diameter;
	}

	void set_lai(double _l){
		lai = _l;
	}

	void set_size(double _x, PlantTraits &traits){
		diameter = _x;
		height = traits.hmat * (1 - exp(-geom.a*diameter/traits.hmat));
		crown_area = geom.pic_4a * height * diameter;
		//leaf_area = crown_area * lai;
		sapwood_fraction = height / (diameter * geom.a);	
	}

	double set_state(std::vector<double>::iterator S, PlantTraits &traits){
		lai = *S++;             // must be set first as it is used bt set_size()
		set_size(*S++, traits);
		litter_pool = *S++;
	}


	double dsize_dmass(PlantTraits &traits) const {
		double dh_dd = geom.a * exp(-geom.a*diameter/traits.hmat);
		double dmleaf_dd = traits.lma * lai * geom.pic_4a * (height + diameter*dh_dd);	// LAI variation is accounted for in biomass production rate
		double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
		double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
		double dmroot_dd = (traits.zeta/traits.lma) * dmleaf_dd;

		double dmass_dd = dmleaf_dd + dmtrunk_dd + dmbranches_dd + dmroot_dd;
		return 1/dmass_dd;
	}

	//double dlai_dt(PlantTraits &traits){
		//return 0.05;
	//}

	double dmass_dt_lai(double dL_dt, PlantTraits &traits){
		return dL_dt * crown_area * (traits.lma + traits.zeta);
	}
	double dlai_dt(double _dmass_dt_lai, PlantTraits &traits){
		return _dmass_dt_lai / (crown_area * (traits.lma + traits.zeta));
	}

	double leaf_mass(PlantTraits &traits){
		return crown_area * lai * traits.lma;	
	}

	double root_mass(PlantTraits &traits){
		return crown_area * lai * traits.zeta;	
	}
	
	//double sapwood_mass(PlantTraits &traits){
		//return traits.wood_density*(hvlc/geom.c)*crown_area*geom.eta_l*height;
	//}
		
	double sapwood_mass(PlantTraits &traits){
		return stem_mass(traits)*sapwood_fraction;
	}

	double sapwood_mass_real(PlantTraits &traits){
		return sapwood_mass(traits)*functional_xylem_fraction;
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

	// projected crown area including gaps, for PPA
	// = r(z)^2 = r0^2 q(z)^2 = Ac/qm^2 q(z)^2 
	double crown_area_extent_projected(double z, PlantTraits &traits){
		if (z >= zm()){
			double fq = q(z)/geom.qm;
			return crown_area * fq*fq;
		} 
		else{
			return crown_area;
		}
	}
		
	double crown_area_above(double z, PlantTraits &traits){
		if (z == 0) return crown_area; // shortcut because z=0 is used often
	
		double fq = q(z)/geom.qm;
		if (z >= zm()){
			return crown_area * fq*fq * (1-geom.fg);
		} 
		else{
			return crown_area * (1 - fq*fq * geom.fg);
		}
	}

	//double leaf_area_above(double z, PlantTraits &traits){
		//return crown_area_above(z, traits) * lai;
	//}

	//double calc_optimal_lai(double P0, double E0, double viscosity_water, PlantParameters &par, PlantTraits &traits){
		//double c1 = par.lambda1;
		//double c2 = par.lambda2 * viscosity_water * E0 * geom.c / (traits.K_xylem);
		//std::cout << "c2 = " << c2 << "\n";

		//double k = par.k_light;
		//double a = (1+c1/c2)*exp(1+k*P0/c2);
		
		//return P0/c2 + 1/k - 1/k * lambertw0(a);
	//}

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
			lai = S[5];
			set_size(S[1], traits);

			double dh_dd = geom.a*exp(-geom.a*diameter/traits.hmat);

			//double dL_dt = -0.05*lai; /[>lai;

			double dB_dt = A*crown_area*lai;	// total biomass production rate
			double dL_dt = 0.05; //dlai_dt(traits);
			double dLA_dt = dmass_dt_lai(dL_dt, traits);  // biomass going into leaf area increment
			double dLit_dt = std::max(-dLA_dt, 0.0);  // biomass going into litter (through leaf loss)
			double dG_dt = dB_dt - std::max(dLA_dt, 0.0); // biomass going into geometric growth
			double dD_dt = dsize_dmass(traits) * dG_dt;	// size (diameter) growth rate

			dSdt[0] = dB_dt;	// biomass that goes into allometric increments
			dSdt[1] = dD_dt; 
			dSdt[2] = 1/(geom.a*diameter*diameter)*(diameter*dh_dd - height)*dD_dt;
			
			double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
			double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
			
			double dsap_trunk_dd = traits.wood_density * M_PI / (4*geom.a) * geom.eta_c * (2*diameter*dh_dd + height) * height;
			double dsap_branch_dd = traits.wood_density * M_PI / (8*geom.a) * sqrt(geom.c/geom.a) * (diameter*dh_dd + height) * sqrt(diameter*height);

			dSdt[3] = (S[3] < sapwood_mass(traits))? (dmtrunk_dd + dmbranches_dd) * dD_dt : (dsap_trunk_dd + dsap_branch_dd) * dD_dt;
			dSdt[4] = (dmtrunk_dd + dmbranches_dd - dsap_trunk_dd - dsap_branch_dd) * dD_dt;
			dSdt[5] = dL_dt;
			dSdt[6] = dLit_dt; //(dL_dt < 0)? dLex_dt : 0;

			k_sap = dSdt[3]/sapwood_mass(traits)*dSdt[1];
		};

		std::vector<double> S = {prod, get_size(), sap_frac_ode, sapwood_mass_ode, heart_mass_ode, lai, litter_pool};
		RK4(t, dt, S, derivs);
		//Euler(t, dt, S, derivs);
		litter_pool = S[6];
		heart_mass_ode = S[4];
		sapwood_mass_ode = S[3];
		sap_frac_ode = S[2];
		lai = S[5];
		set_size(S[1], traits);
		prod = S[0];
		functional_xylem_fraction = S[3]/sapwood_mass(traits);
	}


};


} // namespace plant

#endif


