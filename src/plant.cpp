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
	int i;
   	i = par.initFromFile(file);
	i = traits.initFromFile(file);

	coordinateTraits();

	//geometry.init(par, traits);
}


int Plant::coordinateTraits(){
	traits.ll = 1/(0.0286*pow(traits.lma, -1.71));  // Leaf Economics Spectrum (Relationship from Wright et al. 2004)
	traits.p50_leaf = traits.p50_xylem/3.01;        // P50 = Pg88/3 = P50X/3
	
	traits.K_leaf = exp(1.71-8.628*traits.lma)*1e-16;
	
	par.c = 4*par.a*exp(3.957265-0.040063*traits.hmat)/M_PI;
	
	geometry.init(par, traits);
	
	return 0;
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
	std::cout << "  a = " << geometry.geom.a << "\n";
	std::cout << "  c = " << geometry.geom.c << "\n";
	std::cout << "  lma = " << traits.lma << "\n";
	std::cout << "  hmat = " << traits.hmat << "\n";
	std::cout << "  K_leaf = " << traits.K_leaf << "\n";

}

	
}	// namespace plant



