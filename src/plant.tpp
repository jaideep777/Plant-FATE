namespace plant{


// LAI model
template<class Env>
double Plant::dlai_dt(PlantAssimilationResult& res, Env &env){
	double lai_curr = geometry->lai;
	geometry->set_lai(lai_curr + par.dl);
	auto res_plus = assimilator->net_production(env, geometry, par, traits);
	geometry->set_lai(lai_curr);
	
	double dnpp_dL = (res_plus.npp - res.npp)/geometry->crown_area/par.dl;
	double dE_dL = (res_plus.trans - res.trans)/geometry->crown_area/par.dl;

	double dL_dt = par.response_intensity*(dnpp_dL - 0.001*dE_dL);
	
	return dL_dt;
}




// Demographics
template<class Env>
double Plant::p_survival_germination(Env &env){

}


template<class Env>
double Plant::mortality_rate(Env &env){
	return 0.01;
}


template<class Env>
double Plant::fecundity_rate(double mass_dt, Env &env){
	return mass_dt/(4*traits.seed_mass);
}


// Shorthand is used for biomass partitioning into geometric growth and LAI growth
// These two lines are shorthand for the conditional below
// option 1.
// double dmass_dt_nonlai = dmass_dt - std::max(dmass_dt_lai, 0.0);
// double dmass_dt_lit = std::max(-dmass_dt_lai, 0.0);
// option 2.
// double dmass_dt_nonlai, dmass_dt_lit;
// if (dmass_dt_lai > 0){  // if lai is increasing, biomass is partitioned into allometric (geometric) growth and lai growth
//	dmass_dt_nonlai = dmass_dt - dmass_dt_lai;
//	dmass_dt_lit = 0;
// }
// else{                   // if lai is decreasing, biomass is fully allocated to allometric growth, and lost leaves go into litter
//	dmass_dt_nonlai = dmass_dt;
//	dmass_dt_lit = -dmass_dt_lai;
// }


template<class Env>
void Plant::calc_growth_rates(Env &env){
	
	auto res = assimilator->net_production(env, geometry, par, traits);	
	
	// calculate and constrain rate of LAI change
	double max_alloc_lai = par.max_alloc_lai*std::max(res.npp, 0.0); // if npp is negative, there can be no lai increment. if npp is positive, max 10% can be allocated to lai increment
	rates.dlai_dt = dlai_dt(res, env);
	rates.dmass_dt_lai = geometry->dmass_dt_lai(rates.dlai_dt, max_alloc_lai, traits);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?
	
	// total biomass growth rate (geometric + LAI driven)
	rates.dmass_dt_tot = std::max(res.npp, 0.0);  // 0 if npp is negative

	// if lai is increasing, biomass is partitioned into lai growth and remaining components
	// if lai is decreasing, lost biomass goes into litter
	double dmass_dt_nonlai = rates.dmass_dt_tot - std::max(rates.dmass_dt_lai, 0.0);
	rates.dmass_dt_lit = std::max(-rates.dmass_dt_lai, 0.0);

	// fraction of biomass going into reproduction and biomass allocation to reproduction
	double fR = geometry->dreproduction_dmass(par, traits);
	rates.dmass_dt_rep = fR*dmass_dt_nonlai;
	
	//  fraction of biomass going into growth and size growth rate
	double dmass_growth_dmass = 1-fR;
	rates.dmass_dt_growth = dmass_growth_dmass * dmass_dt_nonlai;
	rates.dsize_dt = geometry->dsize_dmass(traits) * rates.dmass_dt_growth; // size growth rate. NOTE: LAI increment is prioritized (already subtracted from npp above)

	// consistency check - see that all biomass allocations add up as expected
	assert(fabs(rates.dmass_dt_tot - (rates.dmass_dt_lai + rates.dmass_dt_lit + rates.dmass_dt_rep + rates.dmass_dt_growth)) < 1e-6);

//	return {rates.dmass_dt_tot, rates.dlai_dt, rates.dsize_dt, rates.dmass_dt_lit, rates.dmass_dt_rep};
}



}	// namespace plant



