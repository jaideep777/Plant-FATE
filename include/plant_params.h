#ifndef PLANT_FATE_PLANT_PARAMS_H_
#define PLANT_FATE_PLANT_PARAMS_H_

#include "initializer.h"
#include <cmath>

namespace plant{

class PlantParameters{
	public:
	// **
	// ** Photosynthesis paramaters  
	// **
	double kphio = 0.087;
	double alpha = 0.1;
	double gamma = 1.96; 

	// **
	// ** Allocation and geometric paramaters  
	// **
	double m, n; // crown shape paramaters
	double a;      // height-diameter allometry
	double c;      // crown area allometry

	//double hv_min;     // minimum huber value

	// **
	// ** Respiration and turnover 
	// **
	double rl;     // leaf respiration rate per unit photosynthetic capacity (r_leaf = rl*vcmax*leaf_area) [kg yr-1]
	double rr;     // fine-root respiration rate [kg yr-1]
	double rs;     // sapwood respiration rate [kg yr-1]

	//double tl;     // leaf lifespan [yr]
	double lr;     // fine root lifespan [yr]
	
	double cbio;
	double y;


	public:
	// precompute some quantities for efficiency
	double eta_c;
	double pic_4a;

	
	public:
	int initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		I.print();
		m = I.getScalar("mc");
		n = I.getScalar("nc");
		a  = I.getScalar("a");
		c  = I.getScalar("c");
		rl  = I.getScalar("rl");
		rr  = I.getScalar("rr");
		rs  = I.getScalar("rs");
		//kl  = I.getScalar("kl");
		lr  = I.getScalar("lr");
		cbio  = I.getScalar("cbio");
		y = I.getScalar("y");

		hv_min = 1/(c*lai_max);
		eta_c = tgamma(mc+1.0)*tgamma(1.0/nc)/(nc*tgamma(mc+1.0+1.0/nc));
		eta_l = tgamma(ml+1.0)*tgamma(1.0/nl)/(nl*tgamma(ml+1.0+1.0/nl));
		pic_4a = M_PI*c/4/a;

		return 0;
	}

	void print(){
		std::cout << "Params:\n";
		std:: cout << "   m = "  << m << "\n";
		std:: cout << "   n = "  << c << "\n";
		std:: cout << "   a  = " << a  << "\n";
		std:: cout << "   c  = " << c  << "\n";
		std:: cout << "   eta_c  = " << eta_c << "\n";
		std:: cout << "   rl  = " << rl  << "\n";
		std:: cout << "   rr  = " << rr  << "\n";
		std:: cout << "   rs  = " << rs  << "\n";
		//std:: cout << "   ll  = " << kl  << "\n";
		std:: cout << "   lr  = " << lr  << "\n";
		std:: cout << "   cbio  = " << cbio  << "\n";
		std:: cout << "   y   = " << y  << "\n";
	}
};

class PlantTraits{
	// fixed (genetic) traits
	public:
	double lma = 0.09;			// leaf mass per leaf area [kg/m2]
	double zeta = 0.14;			// root mass per leaf area [kg/m2]
	double hmat = 20;			// height at maturity [m]
	double seed_mass = 3.8e-5;	// [kg]
	double wood_density = 608;	// [kg/m3]
	
	double p50_leaf = -1.5;
	double K_leaf = 1e-16;
	double b_leaf = 1;

	// variable (plastic) traits
	public:
	//double vcmax = 40*1e-6*86400*365.2524;	// current vcmax [umol CO2 m-2 s-1]
	double lai_opt;		// optimum crown leaf area index
	double lai;			// realized crown LAI 	
	double ll = 2;	// leaf-longevity (as a function of LMA and environment)
};


} // namespace plant

#endif

