#ifndef PLANT_FATE_PLANT_PARAMS_H_
#define PLANT_FATE_PLANT_PARAMS_H_

#include "utils/initializer.h"
#include <cmath>
#include <string>

namespace plant{

class PlantTraits{
	public:
	std::string species_name = "Tectona grandis";
	
	// fixed (genetic) traits
	public:
	double lma;             // leaf mass per leaf area [kg/m2]
	double zeta;            // root mass per leaf area [kg/m2]
	double hmat;            // height at maturity [m]
	double fhmat;           // height at reproductive maturity as fraction of hmat
	double seed_mass;       // [kg]
	double wood_density;    // [kg/m3]
	double p50_xylem;       // Xylem P50 [MPa]
	double K_leaf;          // Leaf conductivity [m]
	double K_xylem;         // Leaf conductivity [m]
	double b_leaf;          // Shape parameter of leaf vulnerabilty curve [-]
	double b_xylem;         // Shape parameter of leaf vulnerabilty curve [-]
	

	// variable (plastic) traits
	public:
	double ll;	// leaf-longevity (as a function of LMA and environment)
	double p50_leaf;		// Leaf hydraulic capacity [MPa]


	public:
	inline int initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();

		lma = I.getScalar("lma");
		zeta = I.getScalar("zeta");
		hmat = I.getScalar("hmat");
		fhmat = I.getScalar("fhmat");
		seed_mass = I.getScalar("seed_mass");
		wood_density = I.getScalar("wood_density");
		p50_xylem = I.getScalar("p50_xylem");
		K_leaf = I.getScalar("K_leaf");
		K_xylem = I.getScalar("K_xylem");
		b_leaf = I.getScalar("b_leaf");	
		b_xylem = I.getScalar("b_xylem");	
	}
	
};


class PlantParameters{
	public:
	// **
	// ** Photosynthesis paramaters  
	// **
	double kphio;
	double alpha;
	double gamma; 

	// **
	// ** Allocation and geometric paramaters  
	// **
	double m, n; // crown shape paramaters
	double a;      // height-diameter allometry
	double c;      // crown area allometry
	double fg;		// upper canopy gap fraction

	// ** LAI optimization
	double response_intensity;	// leaf construction costs
	double max_alloc_lai;       // max fraction of NPP that can be allocated to LAI increment
	double dl;	                // stepsize for profit derivative
	double lai0;                // initial lai

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

	double k_light;		// light extincttion coefficient

	// ** 
	// ** Demographics
	// **
	double a_f1;    // max fractional allocation to reproduction
	double a_f2;    // rate of increase in reproductive investment
	
	double ll_seed;    // longevity of seeds in the seed pool
	
	// **
	// ** Dispersal and germination
	// **
	double Sd;            // probability of survival during dispersal
	double npp_Sghalf;    // required productivity for 0.5 probability of survival during germination
	
	// **
	// ** Mortality
	// **
	double mI; // baseline mortality rate
	double mD, mD_e; // intrinsic diameter-dependent mortality 
	double mS, mS0; // mortality due to carbon starvation
	
	
	public:
	inline int initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		I.print();
		
//		#define GET(x) x = I.getScalar(#_x);
		kphio = I.getScalar("kphio");
		alpha = I.getScalar("alpha");
		gamma = I.getScalar("gamma");
		m = I.getScalar("m");
		n = I.getScalar("n");
		fg = I.getScalar("fg");
		a  = I.getScalar("a");
		c  = I.getScalar("c");
		response_intensity  = I.getScalar("response_intensity");
		max_alloc_lai  = I.getScalar("max_alloc_lai");
		dl  = I.getScalar("lai_deriv_step");
		lai0  = I.getScalar("lai0");
		rl  = I.getScalar("rl");
		rr  = I.getScalar("rr");
		rs  = I.getScalar("rs");
		//kl  = I.getScalar("kl");
		lr  = I.getScalar("lr");
		cbio  = I.getScalar("cbio");
		y = I.getScalar("y");
		k_light = I.getScalar("k_light");
		a_f1 = I.getScalar("a_f1");
		a_f2 = I.getScalar("a_f2");
		ll_seed = I.getScalar("ll_seed")/log(2);

		Sd = I.getScalar("Sd");
		npp_Sghalf = I.getScalar("npp_Sghalf");

		mI = I.getScalar("mI");
		mD = I.getScalar("mD");
		mD_e = I.getScalar("mD_e");
		mS = I.getScalar("mS");
		mS0 = I.getScalar("mS0");

		return 0;
	}

	inline void print(){
		std::cout << "Params:\n";
		std:: cout << "   m = "  << m << "\n";
		std:: cout << "   n = "  << n << "\n";
		std:: cout << "   a  = " << a  << "\n";
		std:: cout << "   c  = " << c  << "\n";
		//std:: cout << "   eta_c  = " << eta_c << "\n";
		std:: cout << "   rl  = " << rl  << "\n";
		std:: cout << "   rr  = " << rr  << "\n";
		std:: cout << "   rs  = " << rs  << "\n";
		//std:: cout << "   ll  = " << kl  << "\n";
		std:: cout << "   lr  = " << lr  << "\n";
		std:: cout << "   cbio  = " << cbio  << "\n";
		std:: cout << "   y   = " << y  << "\n";
	}
};



} // namespace plant

#endif

