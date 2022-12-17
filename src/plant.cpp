#include "plant.h"
using namespace std;

namespace plant{


void Plant::initParamsFromFile(std::string file){
	par.initFromFile(file);
	traits.initFromFile(file);

	//seeds_hist.set_interval(par.T_seed_rain_avg);

	coordinateTraits();

	//geometry.init(par, traits);
}


void Plant::coordinateTraits(){
	traits.ll = 1/(0.0286*pow(traits.lma, -1.71));  // Leaf Economics Spectrum (Relationship from Wright et al. 2004)
	traits.p50_leaf = traits.p50_xylem/3.01;        // P50 = Pg88/3 = P50X/3
	
//	traits.K_leaf = exp(1.71-8.628*traits.lma)*1e-16;
	
	// double c0 = par.c; // default value of c
	// par.c = 4*par.a* exp(log(c0/2000)+3.957265-0.040063*traits.hmat) /M_PI;
	
	//par.n = 1.1+6*(1-pow(0.5,pow(traits.hmat/25,4)));          // 12*(1-exp(-1*traits.hmat/30));

	par.c = exp(8.968 - 2.6397*traits.hmat/50.876);
	par.a = exp(5.886 - 1.4952*traits.hmat/50.876);

	geometry.init(par, traits);
	
}


void Plant::set_size(double x){
	geometry.set_size(x, traits);
}

double Plant::get_biomass() const{
	return geometry.total_mass(traits);
}

void Plant::set_traits(std::vector<double> tvec){
	traits.lma = tvec[0];
	traits.wood_density = tvec[1];
	coordinateTraits();
}

std::vector<double> Plant::get_traits(){
	vector<double> tvec({
		traits.lma,
		traits.wood_density
	});
	return tvec;
}


void Plant::print(){
	cout << "Plant:\n";
	cout << "  height = " << geometry.height << "\n";
	cout << "  diameter = " << geometry.diameter << "\n";
	cout << "  crown_area = " << geometry.crown_area << "\n";
	cout << "  a = " << geometry.geom.a << "\n";
	cout << "  c = " << geometry.geom.c << "\n";
	cout << "  K_leaf = " << traits.K_leaf << "\n";
	cout << "  lma = " << traits.lma << "\n";
	cout << "  hmat = " << traits.hmat << "\n";
	cout << "  wd = " << traits.wood_density << "\n";
	cout << "  p50 = " << traits.p50_xylem << "\n";

}

	
}	// namespace plant




