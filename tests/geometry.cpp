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
	par.print();
	
	G.initGeometry(par.a, par.c, par.m, par.n, par.fg);

	ofstream fout("geometric_growth.txt");
	fout << "i" << "\t"
		 << "height" << "\t"	
		 << "diameter" << "\t"	
		 << "crown_area" << "\t"	
		 << "leaf_area" << "\t"	
		 << "sapwood_fraction" << "\n";	
	for (int i=1; i<100; ++i){
		G.set_size((i/1000.0)*traits.hmat, traits);
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
		 << "assim" << "\t"
		 << "height" << "\t"	
		 << "diameter" << "\t"	
		 << "crown_area" << "\t"	
		 << "leaf_area" << "\t"	
		 << "sapwood_fraction_ode" << "\t"
		 << "sapwood_fraction" << "\t"
		 << "sapwood_mass" << "\t"
		 << "sapwood_turnover_rate" << "\t"
		 << "heartwood_mass" << "\t"
		 << "heartwood_mass_ode" << "\t"
		 << "total_mass" << "\t"
		 << "total_prod" << "\n";
	G.set_size(0.01, traits);
	G.sap_frac_ode = G.sapwood_fraction;
	G.heart_mass_ode = G.heartwood_mass(traits);

	double dt = 0.1; // mass balance approx holds only for dt = 0.0001
	double total_prod = G.total_mass(traits);
	for (double t=0; t<=100; t=t+dt){

		//cout << t << " " << G.total_mass(traits) << " " << total_prod << "\n";
		if (abs(G.total_mass(traits) - total_prod) > 1e-6) return 1;
		if (abs(G.heartwood_mass(traits) - G.heart_mass_ode) > 1e-6) return 1;
		if (abs(G.sapwood_fraction - G.sap_frac_ode) > 1e-6) return 1;
		
		double P = 1;//*max(sin(M_PI/12.0*t), 0.0); //0.75;	// biomass generation rate - kg/m2/yr
		
		fout << t << "\t"
			 << P << "\t"
			 << G.height << "\t"	
			 << G.diameter << "\t"	
			 << G.crown_area << "\t"	
			 << G.leaf_area << "\t"	
			 << G.sap_frac_ode << "\t"	
			 << G.sapwood_fraction << "\t"	
			 << G.sapwood_mass(traits) << "\t"	
			 << G.k_sap << "\t"	
			 << G.heartwood_mass(traits) << "\t"	
			 << G.heart_mass_ode << "\t"	
			 << G.total_mass(traits) << "\t"
			 << total_prod << "\n";
		
		//total_prod += P*G.leaf_area*dt;
		
		G.grow_for_dt(t, dt, total_prod, P, traits);
		//double dhdM = G.dheight_dmass(traits);
		//double h_new = G.height + dhdM*P*G.leaf_area*dt;
		//G.set_height(h_new, traits);

	}
	fout.close();

	cout << "At t = " << 100 << "\n"
		 << "Total biomass    = " << G.total_mass(traits) << "\n"
		 << "Total production = " << total_prod << "\n";
	
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


//# R script to test:
//dat = read.delim("/home/jaideep/codes/plant_fate_ppa/assim.txt")
//dat$heartwood_fraction = 1-dat$sapwood_fraction

//par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

//plot(dat$height~dat$diameter)
//plot(dat$crown_area~I(dat$height*dat$diameter))
//plot(dat$crown_area~dat$height)
//plot(dat$crown_area~dat$diameter)
//plot(dat$leaf_area~dat$crown_area)
//plot(dat$sapwood_fraction~dat$height)
//# plot(sqrt(4*dat$crown_area/pi)~dat$height)

//plot(dat$height~dat$i)
//plot(dat$diameter~dat$i)
//plot(dat$total_mass~dat$i)
//points(dat$total_prod~dat$i, type="l", col="red")




