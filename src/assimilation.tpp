namespace plant{

// **
// ** Gross and Net Assimilation 
// **
template<class _Climate>
phydro::PHydroResult Assimilator::leaf_assimilation_rate(double I0, double fapar, _Climate &clim, PlantParameters &par, PlantTraits &traits){
	phydro::ParCost par_cost(par.alpha, par.gamma);
	phydro::ParPlant par_plant(traits.K_leaf, traits.p50_leaf, traits.b_leaf);
	par_plant.gs_method = phydro::GS_APX;
	auto photo_leaf = phydro::phydro_analytical(clim.tc,       I0,   clim.vpd,  clim.co2,	
												clim.elv,   fapar,  par.kphio,  clim.swp, 
												par.rl, par_plant,   par_cost);
	
	return photo_leaf;	// umol m-2 s-1 
}


template<class Env>
void  Assimilator::calc_plant_assimilation_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//double GPP_plant = 0, Rl_plant = 0, dpsi_avg = 0;
	double fapar = 1-exp(-0.5*G->lai);

	plant_assim.gpp       = 0;
	plant_assim.rleaf     = 0;
	plant_assim.trans     = 0;
	plant_assim.dpsi_avg  = 0;
	plant_assim.vcmax_avg = 0;
	
	double ca_cumm = 0;
//	std::cout << "--- PPA Assim begin ---" << "\n";
	for (int ilayer=0; ilayer <= env.n_layers; ++ilayer){ // for l in 1:layers{	
		double zst = env.z_star[ilayer];
		double I_top = env.clim.ppfd_max * env.canopy_openness[ilayer]; 
		double ca_layer = G->crown_area_above(zst, traits) - ca_cumm;
		
		auto res = leaf_assimilation_rate(I_top, fapar, env.clim, par, traits);
		
		plant_assim.gpp       += (res.a + res.vcmax*par.rl) * ca_layer;
		plant_assim.rleaf     += (res.vcmax*par.rl) * ca_layer;
		plant_assim.trans     += res.e * ca_layer;
		plant_assim.dpsi_avg  += res.dpsi * ca_layer;
		plant_assim.vcmax_avg += res.vcmax * ca_layer;

		ca_cumm += ca_layer;
		
//		std::cout << "h = " << G->height << ", z* = " << zst << ", I = " << env.canopy_openness[ilayer] << ", fapar = " << fapar << ", A = " << (res.a + res.vcmax*par.rl) << " umol/m2/s x " << ca_layer << " m2 = " << (res.a + res.vcmax*par.rl) * ca_layer << "\n"; 
	}
	assert(fabs(ca_cumm/G->crown_area - 1) < 1e-6);
//	std::cout << "CA traversed = " << ca_cumm << " -- " << G->crown_area << "\n";

	// calculate yearly averages in mol/yr	
	double f_light_day = env.clim.ppfd/env.clim.ppfd_max; //0.25; // fraction day that receives max light (x0.5 sunlight hours, x0.5 average over sinusoid)
	double f_growth_yr = 1.0;  // factor to convert daily mean PAR to yearly mean PAR
	double f = f_light_day * f_growth_yr * 86400*365.2524; // s-1 ---> yr-1

	double ca_total = G->crown_area;                   // total crown area
	plant_assim.gpp   *= (f * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/yr --> mol co2/yr --> kg/yr 
	plant_assim.rleaf *= (f * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/yr --> mol co2/yr --> kg/yr 
	plant_assim.trans *= (f * 18e-3);                  // mol h2o/s  ----> mol h2o/yr  --> kg h2o /yr
	plant_assim.dpsi_avg  /= ca_total;                 // MPa
	plant_assim.vcmax_avg /= ca_total;                 // umol/m2/s
	
}


template<class Env>
PlantAssimilationResult Assimilator::net_production(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	plant_assim = PlantAssimilationResult(); // reset plant_assim

	calc_plant_assimilation_rate(env, G, par, traits); // update plant_assim
	   
	plant_assim.rleaf = leaf_respiration_rate(G,par,traits);      // kg yr-1  
	plant_assim.rroot = root_respiration_rate(G, par,traits);     // kg yr-1
	plant_assim.rstem = sapwood_respiration_rate(G, par,traits);  // kg yr-1
	
	plant_assim.tleaf = leaf_turnover_rate(G, par,traits);  // kg yr-1
	plant_assim.troot = root_turnover_rate(G, par,traits);  // kg yr-1
	
	double A = plant_assim.gpp;
	double R = plant_assim.rleaf + plant_assim.rroot + plant_assim.rstem;
	double T = plant_assim.tleaf + plant_assim.troot;

	plant_assim.npp = par.y*(A-R) - T;	// net biomass growth rate (kg yr-1)

//	std::cout << "assim net = " << plant_assim.npp << ", assim_gros	= " << plant_assim.gpp << "\n"; std::cout.flush();
	return plant_assim;
}

} // namespace plant

