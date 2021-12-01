#include "plant.h"

namespace plant{

Plant::Plant(){
	assimilator = new Assimilator();
	geometry = new PlantGeometry();
}

Plant::~Plant(){
	delete assimilator;
	delete geometry;
}


int Plant::initParamsFromFile(std::string file){
	int i = par.initFromFile(file);
	geometry->initGeometry(0.01, par, traits);
}

void Plant::set_size(double x){
	geometry->set_size(x, traits);
}

double Plant::get_biomass(){
	return geometry->total_mass(traits);
}




// demographics
double Plant::p_survival_germination(){

}

double Plant::mortality_rate(){

}

double Plant::fecundity_rate(double mass, PlantTraits &traits){
	return mass/(4*traits.seed_mass);
}



void Plant::print(){
	std::cout << "Plant:\n";
	std::cout << "  height = " << geometry->height << "\n";
	std::cout << "  diameter = " << geometry->diameter << "\n";
	std::cout << "  crown_area = " << geometry->crown_area << "\n";
	std::cout << "  lai = " << geometry->lai << "\n";
}

	
}	// namespace plant


