#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

namespace plant{

class PlantParameters{
	public:
	// **
	// ** Allocation and geometric paramaters  
	// **
	double ml, nl; // vertical leaf distribution paramaters 
	double mc, nc; // crown shape paramaters
	
	double a;      // height-diameter allometry
	double c;      // crown area allometry
	double b;      // bark allometry

	double lai_max;    // maximum crown leaf area index
	double hv_min;     // minimum huber value

	// **
	// ** Respiration and turnover 
	// **
	double rl;     // leaf respiration rate
	double rr;     // fine-root respiration rate
	double rs;     // sapwood respiration rate

	double kl;     // inverse leaf turnover rate (inverse of leaf longevity)
	double kr;     // fine root turnover rate

	public:
	// precompute some quantities for efficiency
	double eta_c;
	double eta_l;
	double pic_4a;

};

class Traits{
	public:
	double lma;	   // leaf mass per leaf area
	double zeta;   // root mass per leaf area
	double hmat;
	double seed_mass;
	double wood_density;
};


} // namespace plant

#endif

