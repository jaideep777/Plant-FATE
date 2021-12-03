#include "plant.h"

namespace plant{

//Plant::Plant(){
//	assimilator = new Assimilator();
//	geometry = new PlantGeometry();
//}

//Plant::~Plant(){
//	delete assimilator;
//	delete geometry;
//}

//Plant::Plant(const Plant& p) {  
//	assimilator = new Assimilator(*p.assimilator);
//	geometry = new PlantGeometry(*p.geometry);
//}


int Plant::initParamsFromFile(std::string file){
	int i = par.initFromFile(file);
	geometry.init(par, traits);
}

void Plant::set_size(double x){
	geometry.set_size(x, traits);
}

double Plant::get_biomass(){
	return geometry.total_mass(traits);
}


void Plant::print(){
	std::cout << "Plant:\n";
	std::cout << "  height = " << geometry.height << "\n";
	std::cout << "  diameter = " << geometry.diameter << "\n";
	std::cout << "  crown_area = " << geometry.crown_area << "\n";
	std::cout << "  lai = " << geometry.lai << "\n";
}

	
}	// namespace plant




