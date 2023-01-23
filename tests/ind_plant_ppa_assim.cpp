#include <iomanip>
#include <fstream>

#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"

using namespace std;

class Environment : public env::Climate, public env::LightEnvironment {
	public:
	void print(double t){
		Climate::print(t);
		LightEnvironment::print();
	}
};

int main(){

	plant::Plant P;
	P.initParamsFromFile("tests/params/p.ini");
	P.traits.n = 1.1;
	P.geometry.lai = 2.5;
	P.geometry.init(P.par, P.traits); // reinit geometry since params changed
	P.set_size(0.4);
		
	P.par.print();
	P.print();

	Environment C;
	C.z_star = {15,10,5,0};
	C.canopy_openness = {1,exp(-0.5*2),exp(-0.5*3.5),exp(-0.5*5)};
	C.n_layers = C.z_star.size()-1;
	
	C.metFile = "tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv";
	C.co2File = "tests/data/CO2_AMB_AmzFACE2000_2100.csv";
	C.init();
	C.print(0);

	auto res = P.assimilator.net_production(C, &P.geometry, P.par, P.traits);	

	return 0;
}


