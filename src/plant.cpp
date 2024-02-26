#include "plant.h"
using namespace std;

namespace plant{

/// TODO: Should init() functions be made constructors, to prevent risk of creating a plant without
/// @brief  This function initializes the plant (traits, par, and geometry) from an Initialzer object
void Plant::init(const PlantParameters &_par, const PlantTraits &_traits){
	par = _par;
	traits = _traits;
	coordinateTraits();
}


/// @brief  This function initializes the plant (traits, par, and geometry) from an Initialzer object
void Plant::init(io::Initializer &I){
	par.init(I);
	traits.init(I);
	init(par, traits);
}

void Plant::initFromFile(std::string file){
	io::Initializer I;
	I.parse(file);
	init(I);
}


void Plant::coordinateTraits(){
	// traits.ll = 1/(0.0286*pow(traits.lma, -1.71));  // Leaf Economics Spectrum (Relationship from Wright et al. 2004)
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
	// traits.lma = *it++;
	traits.wood_density = *it++;
	traits.hmat = *it++;
	// traits.zeta = *it++;
	init(par, traits);
}

std::vector<double> Plant::get_evolvableTraits(){
	vector<double> tvec({
		// traits.lma,
		traits.wood_density,
		traits.hmat
		// traits.zeta
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




