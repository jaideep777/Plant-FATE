namespace plant{

template<class Env>
std::vector<double> Plant::growth_rates(Env &env){
	
	auto res = assimilator->biomass_growth_rate(env, geometry, par, traits);	
	
	double dL_dt = dlai_dt(res, env, par, traits);

	double max_alloc_lai = par.max_alloc_lai*std::max(res.npp, 0.0); // if npp is negative, there can be no lai increment. if npp is positive, max 10% can be allocated to lai increment
	double dmass_dt_lai = geometry->dmass_dt_lai(dL_dt, max_alloc_lai, traits);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?
		
	double dmass_dt = std::max(res.npp, 0.0);  // total mass increment (geometric and lai-driven), 0 if npp is negative

	double dmass_dt_geom = dmass_dt - std::max(dmass_dt_lai, 0.0);	 // biomass change due to allometric growth. if LAI is decreasing, no mass increase due to LAI
	double dlitter_dt = std::max(-dmass_dt_lai, 0.0);	// biomass from leaf loss goes into litter. if LAI is decreasing, leaves lost go into litter
	
	double dmass_growth_dmass = 1-geometry->dreproduction_dmass(par, traits);
	double dsize_dt = geometry->dsize_dmass(traits) * dmass_growth_dmass * dmass_dt_geom; // size growth rate. NOTE: LAI increment is prioritized (already subtracted from npp above)

	return {dmass_dt, dL_dt, dsize_dt, dlitter_dt, (1-dmass_growth_dmass)*dmass_dt_geom};
}


// LAI model
template<class Env>
double Plant::dlai_dt(PlantAssimilationResult& res, Env &env, PlantParameters &par, PlantTraits &traits){
	double lai_curr = geometry->lai;
	geometry->set_lai(lai_curr + par.dl);
	auto res_plus = assimilator->biomass_growth_rate(env, geometry, par, traits);
	geometry->set_lai(lai_curr);
	
	double dnpp_dL = (res_plus.npp - res.npp)/geometry->crown_area/par.dl;
	double dE_dL = (res_plus.trans - res.trans)/geometry->crown_area/par.dl;

	double dL_dt = par.response_intensity*(dnpp_dL - 0.001*dE_dL); //geometry->dlai_dt(traits);
	
	return dL_dt;
}


}	// namespace plant

