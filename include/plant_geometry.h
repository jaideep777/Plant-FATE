#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#include <cmath>

#include "plant_params.h"

namespace plant{

class PlantGeometry{
	public:
	double height;	// height
	double diameter;	// basal diameter
	double crown_area;	// crown area
	double leaf_area;	// leaf area
	double sapwood_fraction = 1;	// sapwood fraction

	//double fl = 1; // current realized fraction of the maximum possible LAI
	double hvlc = 1;

	void set_height(double _h, PlantParameters &par, PlantTraits &traits){
		height = _h;
		double lai = par.lai_max * traits.fl;
		double hv_min = 1/(par.lai_max * par.c);
		double hv = 1/(lai*par.c);
		diameter = -log(1-height/traits.hmat) * traits.hmat/par.a;
		crown_area = par.pic_4a * height * diameter;
		leaf_area = crown_area*lai;
		sapwood_fraction = (hvlc) * height/diameter/par.a;	// FIXME: check LAI variation
	}

	double dheight_dmass(PlantParameters &par, PlantTraits &traits) const {
		double lai = par.lai_max * traits.fl;
		double dd_dh = 1/(par.a*(1-height/traits.hmat));
		double dmleaf_dh = traits.lma*lai*par.pic_4a * (height*dd_dh + diameter);	// FIXME: Carefully check LAI variation
		double dmstem_dh = (par.eta_l*M_PI*traits.wood_density/4) * (2*height*dd_dh + diameter)*diameter;
		double dmroot_dh = (traits.zeta/traits.lma) * dmleaf_dh;

		double dmass_dh = dmleaf_dh + dmstem_dh + dmroot_dh;
		return 1/dmass_dh;
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

};



} // namespace plant

#endif


