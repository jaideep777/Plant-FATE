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

	double hv; // current realized huber value
	double lai;

	void set_height(double _h, PlantParamaters &par, Traits &traits){
		height = _h;
		diameter = -log(1-height/traits.hmat) * traits.hmat/par.a;
		crown_area = par.pic_4a * height * diamater;
		leaf_area = crown_area*lai;
		sapwood_fraction = (hv*lai*par.c) * height/diamater/par.a;
	}

	double dheight_dmass(PlantParamaters &par, Traits &traits) const {
		double dd_dh = 1/(par.a*(1-height/traits.hmat));
		double dmleaf_dh = traits.lma*lai*par.pic_4a * (height*dd_dh + diamater);
		double dmstem_dh = (par.eta_l*M_PI*traits.wood_density/4) * (2*height*dd_dh + diamater)*diamater;
		double dmroot_dh = (traits.zeta/traits.lma) * dmleaf_dh;

		double dmass_dh = dmleaf_dh + dmstem_dh + dmroot_dh;
		return 1/dmass_dh;
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

}



} // namespace plant

#endif


