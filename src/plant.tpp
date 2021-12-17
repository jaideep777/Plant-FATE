namespace plant{


// LAI model
template<class Env>
double Plant::lai_model(PlantAssimilationResult& res, double _dmass_dt_tot, Env &env){
	double lai_curr = geometry.lai;
	geometry.set_lai(lai_curr + par.dl);
	auto res_plus = assimilator.net_production(env, &geometry, par, traits);
	geometry.set_lai(lai_curr);
	
	double dnpp_dL = (res_plus.npp - res.npp)/geometry.crown_area/par.dl;
	double dE_dL = (res_plus.trans - res.trans)/geometry.crown_area/par.dl;

	double dL_dt = par.response_intensity*(dnpp_dL - 0.001*dE_dL);
	
	if (lai_curr < 0.05 || res.npp < 0) dL_dt = 0;  // limit to prevent LAI going negative
	
	// calculate and constrain rate of LAI change
	double max_alloc_lai = par.max_alloc_lai * _dmass_dt_tot; // if npp is negative, there can be no lai increment. if npp is positive, max 10% can be allocated to lai increment
	bp.dmass_dt_lai = geometry.dmass_dt_lai(dL_dt, max_alloc_lai, traits);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?

	return dL_dt;
}



// seed and sapling survival 
template<class Env>
double Plant::p_survival_germination(Env &env){
	auto res = assimilator.net_production(env, &geometry, par, traits);
	double P = std::max(res.npp, 0.0)/geometry.crown_area;
	double P2 = P*P;
	double P2_half = par.npp_Sghalf * par.npp_Sghalf;
//	std::cout << "P_seed = " << P << ", p_germ = " << P2 / (P2 + P2_half) << std::endl;
	return P2 / (P2 + P2_half);
	
	//calc_demographic_rates(env);
	//return exp(-mortality_rate(env)*0.5);  // mortality during germination realized over half a year
}

template<class Env>
double Plant::p_survival_dispersal(Env &env){
	return par.Sd;
}


// Demographics
template<class Env>
double Plant::size_growth_rate(double _dmass_dt_growth, Env &env){
	double dsize_dt = geometry.dsize_dmass(traits) * _dmass_dt_growth; 
	rates.rgr = dsize_dt/geometry.get_size();
	return dsize_dt;
}


template<class Env>
double Plant::mortality_rate(Env &env){
	double D = geometry.diameter;
	double dDs = par.mS0*exp(-rates.rgr*par.mS); //-log(par.mS0 + rates.rgr*par.mS); //exp(-par.mS * bp.dmass_dt_growth/geometry.crown_area); // Falster-like mortality rate
	double dDd = exp(-par.mD_e*log(D)); //0.1/(1+rates.rgr/0.1);
//	std::cout << "H = " << geometry.height << ", RGR = " << rates.rgr << ", Mortality growth-dependent = " << dD2 << "\n";
	//return par.mI + par.mD*dDd*(1+dDs);

//	double wd = (traits.wood_density/1000);
//	double dI = exp(-5);
//	double mu_rgr = exp(-par.mS*rates.rgr);
//	double mu_d   = exp(-0.3*log(D) + 0.1*D - 1.48*wd*wd);
//	return dI*(mu_d + mu_rgr);
	
	double wd = (traits.wood_density/1000);
//	double dI = exp(-5);
//	double mu_rgr = exp(-par.mS*rates.rgr);
//	double mu_d   = exp(-5 - 0.3*log(D) + 0.2*D - 1.48*(wd*wd-0.690*0.690) + par.mS0*exp(-par.mS*bp.dmass_dt_tot));
	double mu_d   = exp(par.c0 + par.clnD*log(D) + par.cD*D + par.cWD*(wd*wd-par.cWD0*par.cWD0) + par.cS0*exp(-par.cS*bp.dmass_dt_tot));
	return mu_d;
	
	//double logit = -5 -1*log(D*1000) - 0.004*D*1000 + -0.3*log(rates.rgr);
	//return 1/(1+exp(-logit));
	 
}


template<class Env>
double Plant::fecundity_rate(double _dmass_dt_rep, Env &env){
	return _dmass_dt_rep/(4*traits.seed_mass);
}


template<class Env>
void Plant::calc_demographic_rates(Env &env){

	res = assimilator.net_production(env, &geometry, par, traits);	
	bp.dmass_dt_tot = std::max(res.npp, 0.0);  // No biomass growth if npp is negative

	// set rates.dlai_dt and bp.dmass_dt_lai
	rates.dlai_dt = lai_model(res, bp.dmass_dt_tot, env);   // also sets rates.dmass_dt_lai
	
	// set all of bp.dmass_dt_xxx
	partition_biomass(bp.dmass_dt_tot, bp.dmass_dt_lai, env); 

	// set core rates
	rates.dsize_dt  = size_growth_rate(bp.dmass_dt_growth, env);
	rates.dmort_dt  = mortality_rate(env); 

	double fec = fecundity_rate(bp.dmass_dt_rep, env);
	rates.dseeds_dt_pool =  -state.seed_pool/par.ll_seed  +  fec * p_survival_dispersal(env);  // seeds that survive dispersal enter seed pool
	rates.dseeds_dt_germ =   state.seed_pool/par.ll_seed;   // seeds that leave seed pool proceed for germincation
	// need to add seed decay
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
void Plant::partition_biomass(double dm_dt_tot, double dm_dt_lai, Env &env){
	
	// NOTE: LAI increment is prioritized (already subtracted from npp above)	// if lai is increasing, biomass is partitioned into lai growth and remaining components
	// if lai is decreasing, lost biomass goes into litter
	double dmass_dt_nonlai = dm_dt_tot - std::max(dm_dt_lai, 0.0);
	bp.dmass_dt_lit = std::max(-dm_dt_lai, 0.0);

	// fraction of biomass going into reproduction and biomass allocation to reproduction
	double fR = geometry.dreproduction_dmass(par, traits);
	bp.dmass_dt_rep = fR*dmass_dt_nonlai;
	
	//  fraction of biomass going into growth and size growth rate
	double dmass_growth_dmass = 1-fR;
	bp.dmass_dt_growth = dmass_growth_dmass * dmass_dt_nonlai;

	// consistency check - see that all biomass allocations add up as expected
	assert(fabs(bp.dmass_dt_tot - (bp.dmass_dt_lai + bp.dmass_dt_lit + bp.dmass_dt_rep + bp.dmass_dt_growth)) < 1e-6);

	//return {rates.dmass_dt_tot, rates.dlai_dt, rates.dsize_dt, rates.dmass_dt_lit, rates.dmass_dt_rep};
}



}	// namespace plant



