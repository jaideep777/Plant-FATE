#include <iostream>
#include <fstream>

#include "plant_params.h"
#include "plant_geometry.h"

using namespace std;

int main(){

	plant::PlantParameters par;
	plant::PlantTraits traits;
	plant::PlantGeometry G;

	par.initFromFile("tests/params/p.ini");

	ofstream fout("geometric_growth.txt");
	fout << "i" << "\t"
		 << "height" << "\t"	
		 << "diameter" << "\t"	
		 << "crown_area" << "\t"	
		 << "leaf_area" << "\t"	
		 << "sapwood_fraction" << "\n";	
	for (int i=1; i<100; ++i){
		G.set_height((i/100.0)*traits.hmat, par, traits);
		fout << i << "\t"
			 << G.height << "\t"	
			 << G.diameter << "\t"	
			 << G.crown_area << "\t"	
			 << G.leaf_area << "\t"	
			 << G.sapwood_fraction << "\n";	
	}
	fout.close();

	fout.open("geometric_growth_2.txt");
	fout << "i" << "\t"
		 << "height" << "\t"	
		 << "diameter" << "\t"	
		 << "crown_area" << "\t"	
		 << "leaf_area" << "\t"	
		 << "sapwood_fraction" << "\t"
		 << "total_mass" << "\t"
		 << "total_prod" << "\n";
	G.set_height(0.1, par, traits);
	double dt = 0.1; // mass balance approx holds only for dt = 0.0001
	double total_prod = G.total_mass(par, traits);
	for (double t=0; t<100; t=t+dt){

		cout << G.total_mass(par, traits) << " " << total_prod << "\n";

		double P = 0.75;	// biomass generation rate - kg/m2/yr
		total_prod += P*G.leaf_area*dt;
		double dhdM = G.dheight_dmass(par, traits);
		double h_new = G.height + dhdM*P*G.leaf_area*dt;
		
		fout << t << "\t"
			 << G.height << "\t"	
			 << G.diameter << "\t"	
			 << G.crown_area << "\t"	
			 << G.leaf_area << "\t"	
			 << G.sapwood_fraction << "\t"	
			 << G.total_mass(par,traits) << "\t"
			 << total_prod << "\n";
		
		G.set_height(h_new, par, traits);

	}
	fout.close();

	return 0;
}

//# R script to test:
//dat = read.delim("/home/jaideep/codes/plant_fate_ppa/geometric_growth.txt")
//dat$heartwood_fraction = 1-dat$sapwood_fraction

//par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
//plot(dat$height~dat$diameter)
//plot(dat$crown_area~I(dat$height*dat$diameter))
//plot(dat$crown_area~dat$height)
//plot(dat$crown_area~dat$diameter)
//plot(dat$leaf_area~dat$crown_area)
//plot(dat$sapwood_fraction~dat$height)
//plot(sqrt(4*dat$crown_area/pi)~dat$height)

//D = dat$diameter
//H = dat$height
//-log(1-H/20)*20/D
