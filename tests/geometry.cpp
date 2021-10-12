#include <iostream>
#include <fstream>

#include "plant_params.h"
#include "plant_geometry.h"

using namespace std;

int main(){

	PlantParamaters par;
	PlantTraits traits;
	PlantGeometry G;

	ofstream fout("geometric_growth.txt");
	fout << "i" << "\t"
		 << "height" << "\t"	
		 << "diameter" << "\t"	
		 << "crown_area" << "\t"	
		 << "leaf_area" << "\t"	
		 << "sapwood_fraction" << "\n";	
	for (int i=0; i<100; ++i){
		G.set_height((i/100.0)*traits.hmat, par, traits);
		fout << i << "\t"
			 << G.height << "\t"	
			 << G.diameter << "\t"	
			 << G.crown_area << "\t"	
			 << G.leaf_area << "\t"	
			 << G.sapwood_fraction << "\n";	
	}

	return 0;
}


