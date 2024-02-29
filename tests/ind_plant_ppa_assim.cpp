#include <iomanip>
#include <fstream>

#include "traits_params.h"
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
		LightEnvironment::print(t);
	}
};

int main(){

	plant::Plant P;
	P.initFromFile("tests/params/p_test_v2.ini");
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

	C.clim_inst.tc = 26.896;         // temperature, deg C
	C.clim_inst.ppfd_max = 1866;   // umol/m2/s
	C.clim_inst.ppfd = 425;        // umol/m2/s
	C.clim_inst.vpd  = 735;        // Pa
	C.clim_inst.co2  = 368.9;      // ppm
	C.clim_inst.elv = 0;           // m.a.s.l
	C.clim_inst.swp = -0.04;       // MPa

	auto res = P.assimilator.net_production(C, &P.geometry, P.par, P.traits);	

	return 0;
}


