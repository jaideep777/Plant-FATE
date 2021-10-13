#ifndef PLANT_FATE_PLANT_PARAMS_H_
#define PLANT_FATE_PLANT_PARAMS_H_

#include "initializer.h"
#include <cmath>

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

	double kl;     // leaf turnover rate (inverse of leaf longevity)
	double kr;     // fine root turnover rate

	public:
	// precompute some quantities for efficiency
	double eta_c;
	double eta_l;
	double pic_4a;

	int initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		I.print();
		ml = I.getScalar("ml");
		nl = I.getScalar("nl");
		mc = I.getScalar("mc");
		nc = I.getScalar("nc");
		a  = I.getScalar("a");
		c  = I.getScalar("c");
		b  = I.getScalar("b");
		lai_max = I.getScalar("lai_max");
		rl  = I.getScalar("rl");
		rr  = I.getScalar("rr");
		rs  = I.getScalar("rs");
		kl  = I.getScalar("kl");
		kr  = I.getScalar("kr");

		hv_min = 1/(c*lai_max);
		eta_c = tgamma(mc+1.0)*tgamma(1.0/nc)/(nc*tgamma(mc+1.0+1.0/nc));
		eta_l = tgamma(ml+1.0)*tgamma(1.0/nl)/(nl*tgamma(ml+1.0+1.0/nl));
		pic_4a = M_PI*c/4/a;
	}

	void print(){
		std::cout << "Params:\n";
		std:: cout << "   ml = " << ml << "\n";
		std:: cout << "   nl = " << nl << "\n";
		std:: cout << "   mc = " << mc << "\n";
		std:: cout << "   nc = " << nc << "\n";
		std:: cout << "   a  = " << a  << "\n";
		std:: cout << "   b  = " << b  << "\n";
		std:: cout << "   c  = " << c  << "\n";
		std:: cout << "   L  = " << lai_max << "\n";
		std:: cout << "   eta_c  = " << eta_c << "\n";
		std:: cout << "   eta_l  = " << eta_l << "\n";
		std:: cout << "   rl  = " << rl  << "\n";
		std:: cout << "   rr  = " << rr  << "\n";
		std:: cout << "   rs  = " << rs  << "\n";
		std:: cout << "   kl  = " << kl  << "\n";
		std:: cout << "   kr  = " << kr  << "\n";
	}
};

class PlantTraits{
	public:
	double lma = 0.09;			// leaf mass per leaf area [kg/m2]
	double zeta = 0.34;			// root mass per leaf area [kg/m2]
	double hmat = 20;			// height at maturity [m]
	double seed_mass = 3.8e-5;	// [kg]
	double wood_density = 608;	// [kg/m3]
};


} // namespace plant

#endif

