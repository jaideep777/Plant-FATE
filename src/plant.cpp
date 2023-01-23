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
	
	traits.c = exp(8.968 - 2.6397*traits.hmat/50.876);
	traits.a = exp(5.886 - 1.4952*traits.hmat/50.876);

	geometry.init(par, traits);
	
}


void Plant::set_size(double x){
	geometry.set_size(x, traits);
}

double Plant::get_biomass() const{
	return geometry.total_mass(traits);
}

void Plant::set_evolvableTraits(std::vector<double> tvec){
	vector<double>::iterator it = tvec.begin();
	traits.lma = *it++;
	traits.wood_density = *it++;
	coordinateTraits();
}

std::vector<double> Plant::get_evolvableTraits(){
	vector<double> tvec({
		traits.lma,
		traits.wood_density
	});
	return tvec;
}


void Plant::print(){
	cout << "Plant size:\n";
	cout << "  height = " << geometry.height << "\n";
	cout << "  diameter = " << geometry.diameter << "\n";
	cout << "  crown_area = " << geometry.crown_area << "\n";
	traits.print();
}

	
}	// namespace plant




