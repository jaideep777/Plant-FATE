// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_


//#include <plant/control.h>
//#include <plant/disturbance.h>
//#include <plant/interpolator.h>
//#include <plant/adaptive_interpolator.h>
//#include <plant/util.h>

/*#include "gsl_interpolator.h"*/
#include "pn_adaptive_interpolator.h"

// JAI: For now, I am creating a simple environment which always returns the same value for canopy openness.

namespace plant {

//class Patch;

class Environment {
	public:
	
//	Patch * patch;
	//double time = 0;
	double fixed_canopy_openness = 1;

	// Disturbance regime paramters
	double shape = 2;
	double mean_interval = 30;
	double scale = 0;
	double p0 = 0;

	AdaptiveInterpolator light_profile;

	public:
	Environment(double openness);
	
	double canopy_openness(double z) const;

	//double patch_survival() const;	// TODO: Rename to something more intuitive

	double patch_survival(double t) const;	// TODO: Rename to something more intuitive
	double patch_age_density(double t) const; 

//  template <typename Function>
//  void compute_light_environment(Function f_canopy_openness, double height_max);

//  template <typename Function>
//  void rescale_light_environment(Function f_canopy_openness, double height_max);

//  double patch_survival() const;
//  double patch_survival_conditional(double time_at_birth) const;

//  void clear();
//  void clear_light_environment();

//  // NOTE: Interface here will change
//  double seed_rain_dt() const;
//  void set_seed_rain_index(size_t x);

//  // * R interface
//  void r_set_seed_rain_index(util::index x);

//  double time;
//  Disturbance disturbance_regime;
//  interpolator::Interpolator light_environment;

	private:
//  std::vector<double> seed_rain;
//  size_t seed_rain_index;
//  interpolator::AdaptiveInterpolator light_environment_generator;
};

//template <typename Function>
//void Environment::compute_light_environment(Function f_canopy_openness,
//                                            double height_max) {
//  light_environment =
//    light_environment_generator.construct(f_canopy_openness, 0, height_max);
//}

//template <typename Function>
//void Environment::rescale_light_environment(Function f_canopy_openness,
//                                            double height_max) {
//  std::vector<double> h = light_environment.get_x();
//  const double min = light_environment.min(), // 0.0?
//    height_max_old = light_environment.max();

//  util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
//  h.back() = height_max; // Avoid round-off error.

//  light_environment.clear();
//  for (auto hi : h) {
//    light_environment.add_point(hi, f_canopy_openness(hi));
//  }
//  light_environment.initialise();
//}

//inline interpolator::AdaptiveInterpolator
//make_interpolator(const Control& control) {
//  using namespace interpolator;
//  return AdaptiveInterpolator(control.environment_light_tol,
//                              control.environment_light_tol,
//                              control.environment_light_nbase,
//                              control.environment_light_max_depth);
//}

//template <typename T>
//Environment make_environment(T p) {
//  return Environment(p.disturbance_mean_interval, p.seed_rain, p.control);
//}

}

#endif
