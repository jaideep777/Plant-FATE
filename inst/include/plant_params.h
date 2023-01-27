#ifndef PLANT_FATE_PLANT_PARAMS_H_
#define PLANT_FATE_PLANT_PARAMS_H_

#include "utils/initializer.h"
#include <io_utils.h>
#include <cmath>
#include <string>
#include <io_utils.h>
#include <cassert>


namespace plant{

/// \ingroup physiology
class PlantTraits{
	public:
	std::string species_name = "Tectona grandis";
	
	// fixed (genetic) traits
	public:
	double lma;             ///< leaf mass per leaf area [kg/m2]
	double zeta;            ///< root mass per leaf area [kg/m2]
	double fcr;             ///< coarse root mass per unit stem mass [-] 
	double hmat;            ///< height at maturity [m]
	double fhmat;           ///< height at reproductive maturity as fraction of hmat
	double seed_mass;       ///< [kg]
	double wood_density;    ///< [kg/m3]
	double p50_xylem;       ///< Xylem hydraulic vulnerability [MPa]
	double K_leaf;          ///< Leaf conductivity [m]
	double K_xylem;         ///< Leaf conductivity [m]
	double b_leaf;          ///< Shape parameter of leaf vulnerabilty curve [-]
	double b_xylem;         ///< Shape parameter of leaf vulnerabilty curve [-]
	double m;               ///< Crown shape flatness at the top
	double n;               ///< Crown top-heaviness 
	double a;               ///< Initial height to diameter ratio 
	double c;               ///< Crown area to sapwood area ratio
	

	// traits set via coordination
	public:
	// double ll;              ///< leaf-longevity (as a function of LMA and environment)
	double p50_leaf;        ///< Leaf hydraulic vulnerability [MPa] (calculated from Xylem P50)


	public:
	inline void init(io::Initializer &I){
		lma = I.getScalar("lma");
		zeta = I.getScalar("zeta");
		fcr = I.getScalar("fcr");
		hmat = I.getScalar("hmat");
		fhmat = I.getScalar("fhmat");
		seed_mass = I.getScalar("seed_mass");
		wood_density = I.getScalar("wood_density");
		p50_xylem = I.getScalar("p50_xylem");
		K_leaf = I.getScalar("K_leaf");
		K_xylem = I.getScalar("K_xylem");
		b_leaf = I.getScalar("b_leaf");	
		b_xylem = I.getScalar("b_xylem");	
		m = I.getScalar("m");
		n = I.getScalar("n");
		a = I.getScalar("a");	
		c = I.getScalar("c");	
	}

	inline void initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		init(I);
	}

	// Just for debugging purposes - to check if 2 plants have the same traits
	inline bool operator == (const PlantTraits& rhs) const {
		return 
		  (this->lma	== rhs.lma &&
		   this->zeta	== rhs.zeta &&
		   this->fcr	== rhs.fcr &&
		   this->hmat	== rhs.hmat &&
		   this->fhmat	== rhs.fhmat &&
		   this->seed_mass	== rhs.seed_mass &&
		   this->wood_density	== rhs.wood_density &&
		   this->p50_xylem	== rhs.p50_xylem &&
		   this->K_leaf	== rhs.K_leaf &&
		   this->K_xylem  == rhs.K_xylem &&
		   this->b_leaf	== rhs.b_leaf &&
		   this->b_xylem == rhs.b_xylem &&
		   this->m	== rhs.m &&
		   this->n	== rhs.n &&
		   this->a	== rhs.a &&
		   this->c	== rhs.c);
	}


	// Changelog:
	// v2: m,n,a,c move to traits from parameters
	inline void save(std::ostream &fout){
		fout << "Traits::v2 ";
		fout << std::quoted(species_name) << ' ';
		fout << std::make_tuple(
					  lma
					, zeta        
					, fcr         
					, hmat        
					, fhmat       
					, seed_mass   
					, wood_density
					, p50_xylem   
					, K_leaf      
					, K_xylem     
					, b_leaf      
					, b_xylem
					, m
					, n
					, a
					, c    
					);
		fout << '\n';
	}


	inline void restore(std::istream &fin){
		std::string s; fin >> s; // discard version number
		assert(s == "Traits::v2");

		fin >> std::quoted(species_name);
		fin >> lma
			>> zeta        
			>> fcr         
			>> hmat        
			>> fhmat       
			>> seed_mass   
			>> wood_density
			>> p50_xylem   
			>> K_leaf      
			>> K_xylem     
			>> b_leaf      
			>> b_xylem
			>> m
			>> n
			>> a
			>> c;
	}

	inline void print(){
		std::cout << "Traits:\n";
		std::cout << "   lma          = " << lma          << '\n';
		std::cout << "   zeta         = " << zeta         << '\n';
		std::cout << "   fcr          = " << fcr          << '\n';
		std::cout << "   hmat         = " << hmat         << '\n';
		std::cout << "   fhmat        = " << fhmat        << '\n';
		std::cout << "   seed_mass    = " << seed_mass    << '\n';
		std::cout << "   wood_density = " << wood_density << '\n';
		std::cout << "   p50_xylem    = " << p50_xylem    << '\n';
		std::cout << "   K_leaf       = " << K_leaf       << '\n';
		std::cout << "   K_xylem      = " << K_xylem      << '\n';
		std::cout << "   b_leaf       = " << b_leaf       << '\n';
		std::cout << "   b_xylem      = " << b_xylem      << '\n';
		std::cout << "   m            = " << m            << '\n';
		std::cout << "   n            = " << n            << '\n';
		std::cout << "   a            = " << a            << '\n';
		std::cout << "   c            = " << c            << '\n';
	}


};


/// \ingroup physiology
class PlantParameters{
	public:
	// **
	// ** Photosynthesis paramaters  
	// **
	double kphio;           ///< Quantum use efficiency
	double alpha;           ///< Cost of maintaining photosynthetic capacity
	double gamma;           ///< Cost of hydraulic risks

	// **
	// ** Allocation and geometric paramaters  
	// **
	double fg;		        ///< upper canopy gap fraction

	// ** LAI optimization
	double Cc;                  ///< leaf construction costs
	double Chyd;                ///< hydraulic costs
	double response_intensity;	///< speed of response to environment
	double max_alloc_lai;       ///< max fraction of NPP that can be allocated to LAI increment
	double dl;	                ///< stepsize for profit derivative
	double lai0;                ///< initial lai
	bool   optimize_lai;

	// ** Leaf Economics
	double les_u;           ///< LES u [dimensionless]
	double les_cc;          ///< Cost of leaf construction and maintenance  [dimensionless]
	double les_k1;          ///< Conversion factor: g biomass / mol CO2 (see cbio below)
	double les_k2;          ///< Conversion factor: (mol-CO2/day) / (umol-CO2/s)
	double les_hT_dH;       ///< Arrhenius temperature response of  [J mol-1]
	double les_hT_c;        ///<   # - 
	double les_molar_R;     ///< Universal gas constant [J mol-1 K-1]

	// **
	// ** Respiration and turnover 
	// **
	double rd;              ///< leaf dark respiration rate per unit photosynthetic capacity (r_leaf = rl*vcmax*leaf_area) [kg yr-1]
	double rr;              ///< fine-root respiration rate [kg yr-1]
	double rs;              ///< sapwood respiration rate [kg yr-1]

	//double lr;              ///< fine root lifespan [yr]

	double cbio;            ///< Biomass expansion factor: kg biomass per mol CO2 
	double y;               ///< Growth respiration factor [-]

	double k_light;		    ///< light extincttion coefficient

	// ** 
	// ** Demographics
	// **
	double a_f1;            ///< max fractional allocation to reproduction
	double a_f2;            ///< rate of increase in reproductive investment
	
	double ll_seed;         ///< longevity of seeds in the seed pool
	
	// **
	// ** Dispersal and germination
	// **
	double Sd;              ///< probability of survival during dispersal
	double npp_Sghalf;      ///< required productivity for 0.5 probability of survival during germination
	
	// **
	// ** Mortality
	// **
	// double mI; // baseline mortality rate
	// double mD, mD_e; // intrinsic diameter-dependent mortality 
	// double mS, mS0; // mortality due to carbon starvation
	
	// double c0, clnD, cD;  // diameter related mortality params
	// double cL, cG;        // light and growth related mortality
	 
	// double cWD, cWD0;     // wood density related mortality params
	// double cS, cS0;       // light related mortality params
	
	double cD0, cD1;
	double m_alpha, m_beta, m_gamma;

	// **
	// ** Patch structure and successsion
	// **
	// double T_seed_rain_avg;
	
	
	public:
	inline void init(io::Initializer &I){
//		#define GET(x) x = I.getScalar(#_x);
		kphio              = I.getScalar("kphio");
		alpha              = I.getScalar("alpha");
		gamma              = I.getScalar("gamma");
		fg                 = I.getScalar("fg");

		Cc                 = I.getScalar("Cc");
		Chyd               = I.getScalar("Chyd");
		response_intensity = I.getScalar("response_intensity");
		max_alloc_lai      = I.getScalar("max_alloc_lai");
		dl                 = I.getScalar("lai_deriv_step");
		lai0               = I.getScalar("lai0");
		optimize_lai       = (I.getScalar("optimize_lai") == 1) ? true:false;

		les_u              = I.getScalar("les_u");
		les_cc             = I.getScalar("les_cc");
		les_k1             = I.getScalar("les_k1");
		les_k2             = I.getScalar("les_k2"); 
		les_hT_dH          = I.getScalar("les_hT_dH");
		les_molar_R        = I.getScalar("les_molar_R");
		les_hT_c           = I.getScalar("les_hT_c");

		rd                 = I.getScalar("rd");
		rr                 = I.getScalar("rr");
		rs                 = I.getScalar("rs");

		cbio               = I.getScalar("cbio");
		y                  = I.getScalar("y");
		k_light            = I.getScalar("k_light");
		a_f1               = I.getScalar("a_f1");
		a_f2               = I.getScalar("a_f2");

		Sd                 = I.getScalar("Sd");
		npp_Sghalf         = I.getScalar("npp_Sghalf");

		cD0                = I.getScalar("cD0");
		cD1                = I.getScalar("cD1");
		m_alpha            = I.getScalar("m_alpha");
		m_beta             = I.getScalar("m_beta");
		m_gamma            = I.getScalar("m_gamma");

	}
	
	
	inline void initFromFile(std::string fname){
		io::Initializer I(fname);
		I.readFile();
		init(I);
	}

	inline void print(){
		std::cout << "Params:\n";
		std:: cout << "   rd    = " << rd  << "\n";
		std:: cout << "   rr    = " << rr  << "\n";
		std:: cout << "   rs    = " << rs  << "\n";
		std:: cout << "   cbio  = " << cbio  << "\n";
		std:: cout << "   y     = " << y  << "\n";
	}

	inline void save(std::ostream &fout){
		fout << "Params::v1 ";
		fout << std::make_tuple(
			  kphio
			, alpha
			, gamma
			, fg
			, Cc
			, Chyd
			, response_intensity
			, max_alloc_lai
			, dl
			, lai0
			, optimize_lai
			, les_u
			, les_cc
			, les_k1
			, les_k2
			, les_hT_dH
			, les_molar_R
			, les_hT_c
			, rd
			, rr
			, rs
			, cbio
			, y
			, k_light
			, a_f1
			, a_f2
			, ll_seed
			, Sd
			, npp_Sghalf
			, cD0
			, cD1
			, m_alpha
			, m_beta
			, m_gamma
				);
		fout << '\n';
	}


	inline void restore(std::istream &fin){
		std::string s; fin >> s; // discard version number
		assert(s == "Params::v1");

		fin >> kphio
			>> alpha
			>> gamma
			>> fg
			>> Cc
			>> Chyd
			>> response_intensity
			>> max_alloc_lai
			>> dl
			>> lai0
			>> optimize_lai
			>> les_u
			>> les_cc
			>> les_k1
			>> les_k2
			>> les_hT_dH
			>> les_molar_R
			>> les_hT_c
			>> rd
			>> rr
			>> rs
			>> cbio
			>> y
			>> k_light
			>> a_f1
			>> a_f2
			>> ll_seed
			>> Sd
			>> npp_Sghalf
			>> cD0
			>> cD1
			>> m_alpha
			>> m_beta
			>> m_gamma
		;
	}

};



} // namespace plant

#endif

