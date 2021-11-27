#ifndef PLANT_FATE_PLANT_PARAMS_H_
#define PLANT_FATE_PLANT_PARAMS_H_

#include "initializer.h"
#include "incbeta.h"
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
	double fg;		// upper canopy gap fraction

	// ** LAI optimization
	double lambda1;	// leaf construction costs
	double lambda2;	// hydraulic costs

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

	public:
	// precompute some quantities for efficiency
	// Precomputed Geometric parameters
	//double eta_c;
	//double pic_4a;
	//double zm_H;
	//double qm;

	
	public:
	int initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		I.print();
		m = I.getScalar("m");
		n = I.getScalar("n");
		fg = I.getScalar("fg");
		a  = I.getScalar("a");
		c  = I.getScalar("c");
		lambda1  = I.getScalar("lambda1");
		lambda2  = I.getScalar("lambda2");
		rl  = I.getScalar("rl");
		rr  = I.getScalar("rr");
		rs  = I.getScalar("rs");
		//kl  = I.getScalar("kl");
		lr  = I.getScalar("lr");
		cbio  = I.getScalar("cbio");
		y = I.getScalar("y");
		k_light = I.getScalar("k_light");

		//eta_c = 0.33; //tgamma(m+1.0)*tgamma(1.0/n)/(nc*tgamma(m+1.0+1.0/n));

		//pic_4a = M_PI*c/(4*a);

		//zm_H = pow((n-1)/(m*n-1), 1/n);
		//qm = m*n * pow((n-1)/(m*n-1), 1-1/n) * pow((m-1)*n/(m*n-1), m-1);

		//std::cout << "m = " << m << ", n = " << n << ", zm/H = " << zm_H << ", qm = " << qm << "\n";
		
		//eta_c = zm_H - m*m*n/(qm*qm) * beta(2-1/n, 2*m-1) * (incbeta(2-1/n, 2*m-1, (n-1)/(m*n-1)) - (1-fg)); 
		
		return 0;
	}

	void print(){
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

class PlantTraits{
	// fixed (genetic) traits
	public:
	double lma = 0.09;			// leaf mass per leaf area [kg/m2]
	double zeta = 0.14;			// root mass per leaf area [kg/m2]
	double hmat = 20;			// height at maturity [m]
	double seed_mass = 3.8e-5;	// [kg]
	double wood_density = 608;	// [kg/m3]
	
	double p50_leaf = -1.5;		// Leaf hydraulic capacity [MPa]
	double K_leaf = 1e-16;		// Leaf conductivity [m]
	double K_xylem = 2e-16;		// Leaf conductivity [m]
	double b_leaf = 1;			// Shape parameter of leaf vulnerabilty curve [-]

	// variable (plastic) traits
	public:
	//double lai_opt = 1.8;		// optimum crown leaf area index
	//double lai = 1.8;			// realized crown LAI 	
	double ll = 2;	// leaf-longevity (as a function of LMA and environment)
};


} // namespace plant

#endif

