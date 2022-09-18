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
												par.rd, par_plant,   par_cost);
	
	return photo_leaf;	// umol m-2 s-1 
}


template<class Env>
void  Assimilator::calc_plant_assimilation_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//double GPP_plant = 0, Rl_plant = 0, dpsi_avg = 0;
	double fapar = 1-exp(-par.k_light*G->lai);
	bool by_layer = false;
	
	plant_assim.gpp        = 0;
	plant_assim.rleaf      = 0;
	plant_assim.trans      = 0;
	plant_assim.dpsi_avg   = 0;
	plant_assim.vcmax_avg  = 0;
	plant_assim.vcmax25_avg = 0;
	plant_assim.mc_avg     = 0;
	plant_assim.gs_avg     = 0;
	plant_assim.c_open_avg = 0;
	
	double ca_cumm = 0;
	//std::cout << "--- PPA Assim begin ---" << "\n";
	for (int ilayer=0; ilayer <= env.n_layers; ++ilayer){ // for l in 1:layers{	
		double zst = env.z_star[ilayer];
		double ca_layer = G->crown_area_above(zst, traits) - ca_cumm;
		//std::cout << "h = " << G->height << ", z* = " << zst << ", I = " << env.canopy_openness[ilayer] << ", fapar = " << fapar << /*", A = " << (res.a + res.vcmax*par.rd) << " umol/m2/s x " <<*/ ", ca_layer = " << ca_layer << /*" m2 = " << (res.a + res.vcmax*par.rd) * ca_layer << ", vcmax = " << res.vcmax <<*/ "\n"; 
		
		if (by_layer == true){
			double I_top = env.clim.ppfd_max * env.canopy_openness[ilayer]; 
			auto res = leaf_assimilation_rate(I_top, fapar, env.clim, par, traits);
			plant_assim.gpp        += (res.a + res.vcmax*par.rd) * ca_layer;
			plant_assim.rleaf      += (res.vcmax*par.rd) * ca_layer;
			plant_assim.trans      += res.e * ca_layer;
			plant_assim.dpsi_avg   += res.dpsi * ca_layer;
			plant_assim.vcmax_avg  += res.vcmax * ca_layer;
			plant_assim.gs_avg     += res.gs * ca_layer;
			plant_assim.vcmax25_avg += res.vcmax25 * ca_layer;
			plant_assim.mc_avg     += res.mc * ca_layer;
		}
		
		plant_assim.c_open_avg += env.canopy_openness[ilayer] * ca_layer;
		ca_cumm += ca_layer;
		
	}
	assert(fabs(ca_cumm/G->crown_area - 1) < 1e-6);
	double ca_total = G->crown_area;                   // total crown area
	plant_assim.c_open_avg /= ca_total;                // unitless
	if (by_layer == true){
		plant_assim.dpsi_avg   /= ca_total;                // MPa
		plant_assim.vcmax_avg  /= ca_total;                // umol CO2/m2/s
		plant_assim.gs_avg     /= ca_total;                // mol CO2/m2/s
		plant_assim.vcmax25_avg /= ca_total;               // umol CO2/m2/s
		plant_assim.mc_avg     /= ca_total;                // unitless
		//std::cout << "--- total (by layer) \n";
		//std::cout << "h = " << G->height << ", nz* = " << env.n_layers << ", I = " << plant_assim.c_open_avg << ", fapar = " << fapar << ", A = " << plant_assim.gpp/ca_total << " umol/m2/s x " << ca_total << " = " << plant_assim.gpp << ", vcmax_avg = " << plant_assim.vcmax_avg << "\n"; 
	}

	if (by_layer == false){
		double I_top = env.clim.ppfd_max * plant_assim.c_open_avg;
		auto res = leaf_assimilation_rate(I_top, fapar, env.clim, par, traits);
		plant_assim.gpp        = (res.a + res.vcmax*par.rd) * ca_total;
		plant_assim.rleaf      = (res.vcmax*par.rd) * ca_total;
		plant_assim.trans      = res.e * ca_total;
		plant_assim.dpsi_avg   = res.dpsi;
		plant_assim.vcmax_avg  = res.vcmax;
		plant_assim.gs_avg     = res.gs;
		plant_assim.vcmax25_avg = res.vcmax25;
		plant_assim.mc_avg     = res.mc;

		//std::cout << "--- total (by avg light)\n";
		//std::cout << "h = " << G->height << ", nz* = " << env.n_layers << ", I = " << plant_assim.c_open_avg << ", fapar = " << fapar << ", A = " << plant_assim.gpp/ca_total << " umol/m2/s x " << ca_total << " = " << plant_assim.gpp << ", vcmax_avg = " << plant_assim.vcmax_avg << "\n"; 
	}
	//std::cout << "---\nCA traversed = " << ca_cumm << " -- " << G->crown_area << "\n";

	// calculate yearly averages in mol/yr	
	// the factor 1.18 accounts for the non-linearity in the instantaneous sub-daily response in the P-hydro model
	double f_light_day = 1.18*env.clim.ppfd/env.clim.ppfd_max; //0.25; // fraction day that receives max light (x0.5 sunlight hours, x0.5 average over sinusoid)
	double f_growth_yr = 1.0;  // factor to convert daily mean PAR to yearly mean PAR
	double f = f_light_day * f_growth_yr * 86400*365.2524; // s-1 ---> yr-1

	plant_assim.gpp   *= (f * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/yr --> mol co2/yr --> kg/yr 
	plant_assim.npp   *= (f * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/yr --> mol co2/yr --> kg/yr 
	plant_assim.rleaf *= (f * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/yr --> mol co2/yr --> kg/yr 
	plant_assim.trans *= (f * 18e-3);                  // mol h2o/s  ----> mol h2o/yr  --> kg h2o /yr
	
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

	plant_assim.npp = par.y*(A-R) - T; // net biomass growth rate (kg yr-1)

	// if (G->height > 15) std::cout << "h/A = " << G->height << " / " << A/G->crown_area << std::endl;
	// if (env.n_layers > 1 && G->height < 5) std::cout << "h/L/ml/mr | A/R/T/Vc = " << G->height << " / " << G->lai << " / " << G->leaf_mass(traits) << " / " << G->root_mass(traits) << " | " << A << " / " << R << " / " << T << " / " << plant_assim.vcmax_avg << "\n"; 
	// std::cout.flush();
	return plant_assim;
}

} // namespace plant




