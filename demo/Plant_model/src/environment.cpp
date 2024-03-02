#include "environment.h"

#include <vector>
#include <iostream>

//#include <plant/parameters.h>

namespace plant {


Environment::Environment(double openness){
	fixed_canopy_openness = openness;
	//light_profile.splineType = Spline::CONSTRAINED_CUBIC;

	// set up disturbance regime
	scale = pow(tgamma(1.0/shape)/shape/mean_interval, shape);
	p0 = shape*pow(scale, 1.0 / shape) / tgamma(1.0 / shape);

	// By default, construct an empty environment with constant light profile
	std::vector <double> x = {0,10,50,100,1000};
	std::vector <double> y = {openness, openness, openness, openness, openness};

	auto f = [](double x){return 1;}; 
   	light_profile.construct(f, 0, 1000);	// create a fixed empty environment by default
}


double Environment::canopy_openness(double z) const {
//  const bool within_canopy = height <= light_environment.max();
//  return within_canopy ? light_environment.eval(height) : 1.0;
//	return 1;
	// FIXME: Need a check here (instead of in the interpolator)
	return light_profile.eval(z);  // patch->canopy_openness(z); //  fixed_canopy_openness; //
}

//double Environment::canopy_openness(double height) const {
////  const bool within_canopy = height <= light_environment.max();
////  return within_canopy ? light_environment.eval(height) : 1.0;
//	return fixed_canopy_openness;
//}

// Computes the probability of survival from 0 to time.
//double Environment::patch_survival() const {
	//return exp(-scale * pow(time, shape));
//}

// Computes the probability of survival from 0 to t.
double Environment::patch_survival(double t) const {
	return exp(-scale * pow(t, shape));
}

double Environment::patch_age_density(double t) const {
	return p0 * patch_survival(t);
}

//// Computes the probability of survival from time_at_birth to time, by
//// conditioning survival over [0,time] on survival over
//// [0,time_at_birth].
//double Environment::patch_survival_conditional(double time_at_birth) const {
//  return disturbance_regime.pr_survival_conditional(time, time_at_birth);
//}

//// Reset the environment.
//void Environment::clear() {
//  time = 0.0;
//  clear_light_environment();
//}

//void Environment::clear_light_environment() {
//  light_environment.clear();
//}

//double Environment::seed_rain_dt() const {
//  if (seed_rain.empty()) {
//    Rcpp::stop("Cannot get seed rain for empty environment");
//  }
//  return seed_rain[seed_rain_index];
//}

//void Environment::set_seed_rain_index(size_t x) {
//  seed_rain_index = x;
//}

//// * R interface
//void Environment::r_set_seed_rain_index(util::index x) {
//  set_seed_rain_index(x.check_bounds(seed_rain.size()));
//}

}
