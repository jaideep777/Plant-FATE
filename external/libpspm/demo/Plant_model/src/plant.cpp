#include "plant.h"
//#include "../include/plant_parameters.h"

#include <cmath>
#include <functional>
#include <vector>

using namespace std;

//PlantParameters par;

namespace plant {

pn::Integrator Plant::plantIntegrator;

//PlantParameters par;

void Plant::initParameters(){
	// * Individual allometry
	// Canopy shape parameter (extra calculation here later)
	par.eta    = trait_eta; //12.0; // [dimensionless]
	// Ratio sapwood area area to leaf area
	par.theta  = 1.0/4669; // [dimensionless]
	// Height - leaf mass scaling
	par.a_l1   = 5.44; // height with 1m2 leaf [m]
	par.a_l2   = 0.306; // dimensionless scaling of height with leaf area
	// Root mass per leaf area
	par.a_r1   = 0.07;  //[kg / m]
	// Ratio of bark area : sapwood area
	par.a_b1   = 0.17; // [dimensionless]

	// * Production
	// Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
	// =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
	//    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
	// Hard coded in value of lma here so that this value doesn't change
	// if that trait changes above.
	par.r_l   = 21000*0.00187 / lma; // JAI: Should be 39.27/lma; 
	// Root respiration per mass [mol CO2 / yr / kg]
	par.r_r   = 217.0;
	// Sapwood respiration per stem mass  [mol CO2 / yr / kg]
	// = respiration per volume [mol CO2 / m3 / yr]
	// /  wood density [kg/m3]
	par.r_s   = 4012.0 / rho;
	// Bark respiration per stem mass
	// assumed to be twice rate of sapwood
	// (NOTE that there is a re-parametrisation here relative to the paper
	// -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
	// fixed multiplication by 2)
	par.r_b    = 2.0 * par.r_s;
	// Carbon conversion parameter
	par.a_y    = 0.7;
	// Constant converting assimilated CO2 to dry mass [kg / mol]
	// (12E-3 / 0.49)
	par.a_bio  = 2.45e-2;
	// Leaf turnover [/yr]
	par.k_l    = 0.4565855*pow(lma/0.1978791, -1.71); // 0.4565855;	// JAI: Changes with LAI
	// Bark turnover [/yr]
	par.k_b    = 0.2;
	// Sapwood turnover [/yr]
	par.k_s    = 0.2;
	// Root turnover [/yr]
	par.k_r    = 1.0;
	// Parameters of the hyperbola for annual LRC
	par.a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
	par.a_p2   = 0.204716166503633; // [dimensionless]

	// * Seed production
	// Accessory cost of reproduction
	par.a_f3  = 3.0 *  omega; // [kg per seed]

	// Maximum allocation to reproduction
	par.a_f1  = 1.0; //[dimensionless]
	// Size range across which individuals mature
	par.a_f2  = 50; // [dimensionless]

	// * Mortality parameters
	// Probability of survival during dispersal
	par.S_D      = 0.25; // [dimensionless]
	// Parameter for seedling survival
	par.a_d0     = 0.1; //[kg / yr / m2]
	// Baseline for intrinsic mortality
	par.d_I      = 0.01; // [ / yr]
	par.c_d0     = 0.52;
	par.c_d1     = 6.5e-3;
	// Baseline rate for growth-related mortality
	par.a_dG1    = 5.5; // [ / yr]
	// Risk coefficient for dry mass production (per area)
	par.a_dG2    = 20.0;// [yr m2 / kg ]

	par.eta_c = 1 - 2/(1 + par.eta) + 1/(1 + 2*par.eta);
	// NOTE: Also pre-computing, though less trivial
	par.height_0    = height_seed();
	par.area_leaf_0 = area_leaf(par.height_0);
}



Plant::Plant(){
	initParameters();

	vars.height = par.height_0; //0.3257146; //0.3920458; //0.3441948;
	vars.area_leaf = par.area_leaf_0; 
}


double Plant::height() const {
	return vars.height;
}

double Plant::height_dt() const {
	return vars.height_dt;
}

void Plant::set_height(double x) {
	vars.height    = x;
	vars.area_leaf = area_leaf(x);
}

double Plant::mortality() const {
	return vars.mortality;
}

double Plant::mortality_dt() const {
	return vars.mortality_dt;
}

void Plant::set_mortality(double x) {
	vars.mortality = x;
}

double Plant::fecundity() const {
	return vars.fecundity;
}

double Plant::fecundity_dt() const {
	return vars.fecundity_dt;
}

void   Plant::set_fecundity(double x) {
	vars.fecundity = x;
}

double Plant::area_heartwood() const {
	return vars.area_heartwood;
}

double Plant::area_heartwood_dt() const {
	return vars.area_heartwood_dt;
}

void   Plant::set_area_heartwood(double x) {
	vars.area_heartwood = x;
}

double Plant::mass_heartwood() const {
	return vars.mass_heartwood;
}

double Plant::mass_heartwood_dt() const {
	return vars.mass_heartwood_dt;
}

void   Plant::set_mass_heartwood(double x) {
	vars.mass_heartwood = x;
}




// [eqn 2] area_leaf (inverse of [eqn 3])
double Plant::area_leaf(double height) const {
  return pow(height / par.a_l1, 1.0 / par.a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
double Plant::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double Plant::area_sapwood(double area_leaf) const {
  return area_leaf * par.theta;
}

double Plant::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * par.eta_c * rho;
}

// [eqn 5] area and mass of bark
double Plant::area_bark(double area_leaf) const {
  return par.a_b1 * area_leaf * par.theta;
}

double Plant::mass_bark(double area_bark, double height) const {
  return area_bark * height * par.eta_c * rho;
}

double Plant::area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double Plant::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double Plant::mass_root(double area_leaf) const {
  return par.a_r1 * area_leaf;
}

// [eqn 8] Total mass
double Plant::mass_live(double mass_leaf, double mass_bark, double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double Plant::mass_total(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_heartwood,
                            double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double Plant::mass_above_ground(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root; // FIXME: Seems wrong: --root-- *heartwood
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by par.a_y*par.a_bio.
double Plant::respiration(double mass_leaf, double mass_sapwood,
                             double mass_bark, double mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

double Plant::respiration_leaf(double mass) const {
  return par.r_l * mass;
}

double Plant::respiration_bark(double mass) const {
  return par.r_b * mass;
}

double Plant::respiration_sapwood(double mass) const {
  return par.r_s * mass;
}

double Plant::respiration_root(double mass) const {
  return par.r_r * mass;
}

// [eqn 14] Total turnover
double Plant::turnover(double mass_leaf, double mass_bark,
                          double mass_sapwood, double mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

double Plant::turnover_leaf(double mass) const {
  return par.k_l * mass;
}

double Plant::turnover_bark(double mass) const {
  return par.k_b * mass;
}

double Plant::turnover_sapwood(double mass) const {
  return par.k_s * mass;
}

double Plant::turnover_root(double mass) const {
  return par.k_r * mass;
}

//// [eqn 15] Net production
////
//// NOTE: Translation of variable names from the Falster 2011.  Everything
//// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
//double Plant::net_mass_production_dt_A(double assimilation, double respiration,
//                                double turnover) const {
//  return par.a_bio * par.a_y * (assimilation - respiration) - turnover;
//}

// [eqn 16] Fraction of production allocated to reproduction
double Plant::fraction_allocation_reproduction(double height) const {
  return par.a_f1 / (1.0 + exp(par.a_f2 * (1.0 - height / hmat)));
}

// Fraction of production allocated to growth
double Plant::fraction_allocation_growth(double height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
double Plant::fecundity_dt(double net_mass_production_dt,
                               double fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (omega + par.a_f3);
}

double Plant::darea_leaf_dmass_live(double area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double Plant::dheight_darea_leaf(double area_leaf) const {
  return par.a_l1 * par.a_l2 * pow(area_leaf, par.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double Plant::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double Plant::dmass_sapwood_darea_leaf(double area_leaf) const {
  return rho * par.eta_c * par.a_l1 * par.theta * (par.a_l2 + 1.0) * pow(area_leaf, par.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double Plant::dmass_bark_darea_leaf(double area_leaf) const {
  return par.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double Plant::dmass_root_darea_leaf(double /* area_leaf */) const {
  return par.a_r1;
}

// Growth rate of basal diameter_stem per unit time
double Plant::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double Plant::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * par.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double Plant::area_heartwood_dt(double area_leaf) const {
  return par.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double Plant::area_bark_dt(double area_leaf_dt) const {
  return par.a_b1 * area_leaf_dt * par.theta;
}

// Growth rate of stem basal area per unit time
double Plant::area_stem_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double Plant::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double Plant::mass_root_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

double Plant::mass_live_dt(double fraction_allocation_reproduction,
                               double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double Plant::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double Plant::mass_above_ground_dt(double area_leaf,
                                       double fraction_allocation_reproduction,
                                       double net_mass_production_dt,
                                       double mass_heartwood_dt,
                                       double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

double Plant::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double Plant::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

double Plant::height_given_mass_leaf(double mass_leaf) const {
  return par.a_l1 * pow(mass_leaf / lma, par.a_l2);
}

double Plant::mortality_dt(double productivity_area,
                              double cumulative_mortality) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (!isinf(cumulative_mortality)){ // < 1e20) {//(R_FINITE(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    // return 0.0;		
    // Resolved [FIXME]: JAI: Shoudnt we return some very large number here so that density reduces via the differential equation?
    return mortality_growth_independent_dt()*100; // JAI fix: return 100x the baseline mortality rate
  }
}

double Plant::mortality_growth_independent_dt() const {
  return par.d_I;
  //return par.c_d0*exp(-par.c_d1*rho);
}

double Plant::mortality_growth_dependent_dt(double productivity_area) const {
  return par.a_dG1 * exp(-par.a_dG2 * productivity_area);
}


double Plant::area_leaf_above(double z, double height,
                                 double area_leaf) const {
  return area_leaf * Q(z, height);	
}

// [eqn  9] Probability density of leaf area at height `z`
double Plant::q(double z, double height) const {
	if (z == 0 || z > height) return 0;
	else{
		//const double tmp = pow(z / height, par.eta);
	    const double tmp = exp(log(z / height)*par.eta);
		return 2 * par.eta * (1 - tmp) * tmp / z;
	}
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double Plant::Q(double z, double height) const {
  if (z > height) return 0;
  if (z == 0) return 1;
  
//  const double tmp = 1.0-pow(z / height, par.eta);
  const double tmp = 1.0-exp(log(z / height)*par.eta);
  return tmp * tmp;
}

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
double Plant::Qp(double x, double height) const { // x in [0,1], unchecked.   // TODO: Better to rename Qinv ??
  return pow(1 - sqrt(x), (1/par.eta)) * height;
}



// Prints root of func(x) with error of EPSILON 
template <typename Function>
double bisection(Function func, double a, double b) { 
    if (func(a) * func(b) >= 0){ 
		std::cout << "You have not assumed right a and b\n"; 
        return 0; 
    } 
  
    double c = a; 
    while ((b-a) >= 1e-6) { 
        // Find middle point 
        c = (a+b)/2; 
  
        // Check if middle point is root 
        if (func(c) == 0.0) break; 
        // Decide the side to repeat the steps 
        else if (func(c)*func(a) < 0)   b = c; 
        else a = c; 
    } 
    cout << "The value of root is : " << c << "\n"; 
	return c;
} 


//// The aim is to find a plant height that gives the correct seed mass.
double Plant::height_seed(void) const {

//  // Note, these are not entirely correct bounds. Ideally we would use height
//  // given *total* mass, not leaf mass, but that is difficult to calculate.
//  // Using "height given leaf mass" will expand upper bound, but that's ok
//  // most of time. Only issue is that could break with obscure parameter
//  // values for LMA or height-leaf area scaling. Could instead use some
//  // absolute maximum height for new seedling, e.g. 1m?
	
  const double
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(omega);

  const double tol = 1e-6; //control.plant_seed_tol;
  const size_t max_iterations = 100; //control.plant_seed_iterations;

  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };

  return bisection(target, h0, h1);

//  return util::uniroot(target, h0, h1, tol, max_iterations);
}

std::ostream& operator<<(std::ostream& os, const Plant& p){
	os << "Plant: " << std::endl;
	os << "> Leaf Area = " << p.vars.area_leaf << std::endl;
	os << "> State\tRate" << std::endl;
	os << "ht : " << p.vars.height            << "\t" << p.vars.height_dt         << std::endl;
	os << "mo : " << p.vars.mortality         << "\t" << p.vars.mortality_dt      << std::endl;
	os << "fe : " << p.vars.fecundity         << "\t" << p.vars.fecundity_dt      << std::endl;
	os << "ha : " << p.vars.area_heartwood    << "\t" << p.vars.area_heartwood_dt << std::endl;
	os << "hm : " << p.vars.mass_heartwood    << "\t" << p.vars.mass_heartwood_dt << std::endl;
	return os;	
}


// New ODE interface
void Plant::set_state(const double * full_state_array){
	vars.height         = full_state_array[0];
	vars.mortality      = full_state_array[1];
	vars.fecundity      = full_state_array[2];
	vars.area_heartwood = full_state_array[3];
	vars.mass_heartwood = full_state_array[4];
	
	vars.area_leaf = area_leaf(vars.height);
}

vector <double> Plant::get_state(){
	vector <double> state = {vars.height, 
							 vars.mortality, 
							 vars.fecundity, 
							 vars.area_heartwood, 
							 vars.mass_heartwood}; 
	return state;
}

void Plant::get_rates(double * all_rates_array){
	all_rates_array[0] = vars.height_dt;
	all_rates_array[1] = vars.mortality_dt;
	all_rates_array[2] = vars.fecundity_dt;
	all_rates_array[3] = vars.area_heartwood_dt;
	all_rates_array[4] = vars.mass_heartwood_dt;
}


// [Appendix S6] Per-leaf photosynthetic rate.
// Here, `x` is openness, ranging from 0 to 1.
double Plant::assimilation_leaf(double x) const {
  return par.a_p1 * x / (x + par.a_p2);
}


//void Plant::step_by(double Dt, Environment &env){

//	double t0 = 0, tf = Dt, dt = 0.05;
//	size_t nsteps = (tf-t0)/dt;

//	for (size_t i=0; i < nsteps; ++i){
//		compute_vars_phys(env, false);
//		
//		// update state
//		vars.height += vars.height_dt * dt;
//		vars.mortality += vars.mortality_dt * dt;
//		vars.fecundity += vars.fecundity_dt * dt;
//		vars.area_heartwood += vars.area_heartwood_dt * dt;
//		vars.mass_heartwood += vars.mass_heartwood_dt * dt;

//		// update other variables
//		vars.area_leaf = area_leaf(vars.height);
//		
//	}

//}



}
