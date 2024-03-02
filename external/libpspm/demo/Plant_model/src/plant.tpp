

// [eqn 12] Gross annual CO2 assimilation
//
// NOTE: In contrast with Daniel's implementation (but following
// Falster 2012), we do not normalise by par.a_y*par.a_bio here.
template <class Environment>
double Plant::assimilation(const Environment& environment,
                                    double height,
                                    double area_leaf,
                                    bool reuse_intervals) {
//  const bool over_distribution = control.plant_assimilation_over_distribution;
  const double x_min = 0.0, x_max = /*over_distribution ? 1.0 :*/ height;

  double A = 0.0;

//  std::function<double(double)> f;
//  if (over_distribution) {
//    f = [&] (double x) -> double {
//      return compute_assimilation_p(x, height, environment);
//    };
//  } else {
    auto f = [&height, &environment, this] (double z) -> double {
      // return compute_assimilation_h(x, height, environment);
	  return assimilation_leaf(environment.canopy_openness(z)) * q(z, height); // JAI: This was assimilation_h
	};
//  }

//  if (control.plant_assimilation_adaptive && reuse_intervals) {
//    A = control.integrator.integrate_with_last_intervals(f, x_min, x_max);
//  } else {
//    A = control.integrator.integrate(f, x_min, x_max);

//		steady_clock::time_point t1 = steady_clock::now();
	  A = plantIntegrator.integrate(f, x_min, x_max);
	  //A = assimilation_leaf(environment.canopy_openness(height)); // JAI: Big leaf model
//		steady_clock::time_point t2 = steady_clock::now();
//	cout << "plant::assimilation [" << duration_cast<duration<double>>(t2 - t1).count() << " sec]" << endl;
		
		
//  }
	
	
	
	//  std::cout << "height = " << height << ", light = " << environment.canopy_openness(height) << ", assimmilation = " << A << " " << area_leaf << std::endl;
  
  return area_leaf * A;
}

// This is used in the calculation of assimilation by
// `compute_assimilation` above; it is the term within the integral in
// [eqn 12]; i.e., A_lf(A_0v, E(z,a)) * q(z,h(m_l))
// where `z` is height.
//double Plant::compute_assimilation_x(double x, double height,
//                                     const Environment& environment) const {
//  if (control.plant_assimilation_over_distribution) {
//    return compute_assimilation_p(x, height, environment);
//  } else {
//    return compute_assimilation_h(x, height, environment);
//  }
//}

//double Plant::compute_assimilation_h(double z, double height,
//                                     const Environment& environment) const {
//  return assimilation_leaf(environment.canopy_openness(z)) * q(z, height);
//}

//double Plant::compute_assimilation_p(double p, double height,
//                                     const Environment& environment) const {
//  return assimilation_leaf(environment.canopy_openness(Qp(p, height)));			// Qp = Qinv
//}

// One shot calculation of net_mass_production_dt
// Used by germination_probability() and scm_vars().
template <class Environment>
double Plant::net_mass_production_dt(const Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals) {
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);

  const double assimilation_ = assimilation(environment, height, area_leaf_, reuse_intervals);
  const double respiration_  = respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_     = turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);

  //return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
  return par.a_bio * par.a_y * (assimilation_ - respiration_) - turnover_;
}



// one-shot update of the scm variables
template <class Environment>
void Plant::compute_vars_phys(const Environment& environment,
                              bool reuse_intervals) {
                     
//  vars.area_leaf = area_leaf(vars.height);	// JAI: Need to update this before calculating rates
  const double net_mass_production_dt_ = net_mass_production_dt(environment, vars.height, vars.area_leaf, reuse_intervals);

  if (net_mass_production_dt_ > 0) {
    const double fraction_allocation_reproduction_ = fraction_allocation_reproduction(vars.height);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(vars.area_leaf);
    const double fraction_allocation_growth_ = fraction_allocation_growth(vars.height);
    const double area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.height_dt = dheight_darea_leaf(vars.area_leaf) * area_leaf_dt;
    vars.fecundity_dt = fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_);

    vars.area_heartwood_dt = area_heartwood_dt(vars.area_leaf);
    const double area_sapwood_ = area_sapwood(vars.area_leaf);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, vars.height);
    vars.mass_heartwood_dt = mass_heartwood_dt(mass_sapwood_);
  } 
  else {
    vars.height_dt         = 0.0;
    vars.fecundity_dt      = 0.0;
    vars.area_heartwood_dt = 0.0;
    vars.mass_heartwood_dt = 0.0;
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.mortality_dt = mortality_dt(net_mass_production_dt_ / vars.area_leaf, vars.mortality);
  //std::cout << "lma | h/g/m/f = " << lma << " | " << vars.height << " -> " << environment.light_profile.eval(vars.height) << " -> " << vars.height_dt << " " << vars.mortality_dt << " " << vars.fecundity_dt << std::endl;
}



// [eqn 20] Survival of seedlings during germination
template <class Environment>
double Plant::germination_probability(const Environment& environment) {
	double height_0 = par.height_0; //height_seed(); //0.344;	// JAI: 0.344 was Temporary provision. // better to pre-compute if height_seed() uses root-finding
  const double net_mass_production_dt_ = net_mass_production_dt(environment, height_0, area_leaf(height_0));
  //std::cout << "net mass = " << height_0 << " " << environment.canopy_openness(0) << " " << net_mass_production_dt_ << std::endl;
  if (net_mass_production_dt_ > 0) {
    const double tmp = par.a_d0 * area_leaf(height_0) / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0);
  } else {
    return 0.0;
  }
}

