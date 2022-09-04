#include "plant_geometry.h"

#include "utils/rk4.h"
#include "utils/incbeta.h"
#include "plant_params.h"

#include <cmath>

namespace plant{

void PlantGeometry::init(PlantParameters &par, PlantTraits &traits){
	geom.m = par.m; geom.n = par.n; 
	geom.a = par.a; geom.c = par.c;
	geom.fg = par.fg;

	geom.pic_4a = M_PI*geom.c/(4*geom.a);

	double m = geom.m, n = geom.n;
	geom.zm_H = pow((n-1)/(m*n-1), 1/n);
	geom.qm = m*n * pow((n-1)/(m*n-1), 1-1/n) * pow((m-1)*n/(m*n-1), m-1);

	geom.eta_c = geom.zm_H - m*m*n/(geom.qm*geom.qm) * beta(2-1/n, 2*m-1) * (incbeta(2-1/n, 2*m-1, (n-1)/(m*n-1)) - (1-geom.fg)); 
	
	// std::cout << "Init Geometry: m = " << m << ", n = " << n << ", zm/H = " << geom.zm_H << ", qm = " << geom.qm << ", eta_c = " << geom.eta_c << "\n";
	
	geom.dmat = -(traits.hmat/geom.a) * log(1-traits.fhmat);

//	lai = par.lai0;
//	set_size(diameter_0, traits);
}

// **
// ** Crown geometry
// **
double PlantGeometry::q(double z){
	if (z > height || z < 0) return 0;
	else{
		double m = geom.m, n = geom.n;
		double zHn_1 = pow(z/height, n-1);
		double zHn   = zHn_1 * z/height;
		return m*n * pow(1 - zHn, m-1) * zHn_1;
	}
}

double PlantGeometry::zm(){
	return geom.zm_H * height;
} 

// projected crown area including gaps, for PPA
// = r(z)^2 = r0^2 q(z)^2 = (Ac/qm^2) q(z)^2 = Ac (q(z)/qm)^2
double PlantGeometry::crown_area_extent_projected(double z, PlantTraits &traits){
	if (z >= zm()){
		double fq = q(z)/geom.qm;
		return crown_area * fq*fq;
	} 
	else{
		return crown_area;
	}
}
	
double PlantGeometry::crown_area_above(double z, PlantTraits &traits){
	if (z == 0) return crown_area; // shortcut because z=0 is used often

	double fq = q(z)/geom.qm;
	if (z >= zm()){
		return crown_area * fq*fq * (1-geom.fg);
	} 
	else{
		return crown_area * (1 - fq*fq * geom.fg);
	}
}


// **
// ** Biomass partitioning
// **
double PlantGeometry::dsize_dmass(PlantTraits &traits) const {
	double dh_dd = geom.a * exp(-geom.a*diameter/traits.hmat);
	double dmleaf_dd = traits.lma * lai * geom.pic_4a * (height + diameter*dh_dd);	// LAI variation is accounted for in biomass production rate
	double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
	double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
	double dmroot_dd = (traits.zeta/traits.lma) * dmleaf_dd;
	double dmcroot_dd = (dmbranches_dd + dmtrunk_dd)*traits.fcr;

	double dmass_dd = dmleaf_dd + dmtrunk_dd + dmbranches_dd + dmroot_dd + dmcroot_dd;
	return 1/dmass_dd;
}

double PlantGeometry::dreproduction_dmass(PlantParameters &par, PlantTraits &traits){
	return par.a_f1 / (1.0 + exp(par.a_f2 * (1.0 - diameter / geom.dmat))); 
}


// **
// ** LAI model
// ** 
double PlantGeometry::dmass_dt_lai(double &dL_dt, double dmass_dt_max, PlantTraits &traits){
	double l2m = crown_area * (traits.lma + traits.zeta); 
	double dm_dt_lai = std::min(dL_dt * l2m, dmass_dt_max);  // biomass change resulting from LAI change  // FIXME: here roots also get shed with LAI. true?
	dL_dt = dm_dt_lai / l2m;
	return dm_dt_lai;
}



// **
// ** Carbon pools
// **
double PlantGeometry::leaf_mass(const PlantTraits &traits) const{
	return crown_area * lai * traits.lma;
}

double PlantGeometry::root_mass(const PlantTraits &traits) const{
	return crown_area * lai * traits.zeta;
}

double PlantGeometry::coarse_root_mass(const PlantTraits &traits) const{
	return stem_mass(traits)*traits.fcr;	
}

//double sapwood_mass(PlantTraits &traits){
	//return traits.wood_density*(hvlc/geom.c)*crown_area*geom.eta_l*height;
//}
	
double PlantGeometry::sapwood_mass(const PlantTraits &traits) const{
	return stem_mass(traits)*sapwood_fraction;
}

double PlantGeometry::sapwood_mass_real(const PlantTraits &traits) const{
	return sapwood_mass(traits)*functional_xylem_fraction;
}

double PlantGeometry::stem_mass(const PlantTraits &traits) const{
	double trunk_mass = traits.wood_density*(M_PI*diameter*diameter/4)*height*geom.eta_c;
	double branch_mass = traits.wood_density * (M_PI*diameter*diameter/12)*height * sqrt((geom.c/geom.a)*(diameter/height));	
	return trunk_mass + branch_mass;
}

double PlantGeometry::heartwood_mass(const PlantTraits &traits) const{
	return stem_mass(traits)*(1-sapwood_fraction);
}

double PlantGeometry::total_mass(const PlantTraits &traits) const{
	return stem_mass(traits)*(1+traits.fcr) + leaf_mass(traits) + root_mass(traits);
}

// **
// ** state manipulations
// **	
double PlantGeometry::get_size() const {
	return diameter;
}

void PlantGeometry::set_lai(double _l){
	lai = _l;
}


void PlantGeometry::set_size(double _x, PlantTraits &traits){
	diameter = _x;
	height = traits.hmat * (1 - exp(-geom.a*diameter/traits.hmat));
	crown_area = geom.pic_4a * height * diameter;
	//leaf_area = crown_area * lai;
	sapwood_fraction = height / (diameter * geom.a);	
}

std::vector<double>::iterator PlantGeometry::set_state(std::vector<double>::iterator S, PlantTraits &traits){
	set_lai(*S++);             // must be set first as it is used bt set_size() - not required any more
	set_size(*S++, traits);
//	litter_pool = *S++;
	return S;
}


// ** 
// ** Simple growth simulator for testing purposes
// ** - simulates growth over dt with constant assimilation rate A
// ** 
void PlantGeometry::grow_for_dt(double t, double dt, double &prod, double &litter_pool, double A, PlantTraits &traits){

	auto derivs = [A, &traits, &litter_pool, this](double t, std::vector<double>&S, std::vector<double>&dSdt){
		set_lai(S[5]);
		set_size(S[1], traits);
		litter_pool = S[6];

		double dh_dd = geom.a*exp(-geom.a*diameter/traits.hmat);

		//double dL_dt = -0.05*lai; /[>lai;

		// ignoring reproduction in these biomass pools
		double dB_dt = A*crown_area;	// total biomass production rate
		double dL_dt = -0.01; //dlai_dt(traits);
		double dLA_dt = dmass_dt_lai(dL_dt, dB_dt, traits);  // biomass going into leaf area increment
		double dLit_dt = std::max(-dLA_dt, 0.0);  // biomass going into litter (through leaf loss)
		double dG_dt = dB_dt - std::max(dLA_dt, 0.0); // biomass going into geometric growth
		double dD_dt = dsize_dmass(traits) * dG_dt;	// size (diameter) growth rate

		dSdt[0] = dB_dt;	// biomass that goes into allometric increments
		dSdt[1] = dD_dt; 
		dSdt[2] = 1/(geom.a*diameter*diameter)*(diameter*dh_dd - height)*dD_dt;
		
		double dmtrunk_dd = (geom.eta_c * M_PI * traits.wood_density / 4) * (2*height + diameter*dh_dd)*diameter;
		double dmbranches_dd = (sqrt(geom.c / geom.a) * M_PI * traits.wood_density / 12) * (2.5*height + 0.5*diameter*dh_dd) * diameter*sqrt(diameter/height); 
		
		double dsap_trunk_dd = traits.wood_density * M_PI / (4*geom.a) * geom.eta_c * (2*diameter*dh_dd + height) * height;
		double dsap_branch_dd = traits.wood_density * M_PI / (8*geom.a) * sqrt(geom.c/geom.a) * (diameter*dh_dd + height) * sqrt(diameter*height);

		dSdt[3] = (S[3] < sapwood_mass(traits))? (dmtrunk_dd + dmbranches_dd) * dD_dt : (dsap_trunk_dd + dsap_branch_dd) * dD_dt;
		dSdt[4] = (dmtrunk_dd + dmbranches_dd - dsap_trunk_dd - dsap_branch_dd) * dD_dt;
		dSdt[5] = dL_dt;
		dSdt[6] = dLit_dt; //(dL_dt < 0)? dLex_dt : 0;

		k_sap = dSdt[3]/sapwood_mass(traits)*dSdt[1];
	};

	std::vector<double> S = {prod, get_size(), sap_frac_ode, sapwood_mass_ode, heart_mass_ode, lai, litter_pool};
	RK4(t, dt, S, derivs);
	//Euler(t, dt, S, derivs);
	litter_pool = S[6];
	heart_mass_ode = S[4];
	sapwood_mass_ode = S[3];
	sap_frac_ode = S[2];
	set_lai(S[5]);
	set_size(S[1], traits);
	prod = S[0];
	functional_xylem_fraction = S[3]/sapwood_mass(traits);
}


} // namespace plant



