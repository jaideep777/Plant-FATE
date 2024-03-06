#include <iomanip>
#include <fstream>

#include "traits_params.h"
#include "plant_architecture.h"
#include "assimilation.h"
#include "plant.h"

#include "climate.h"
#include "light_environment.h"

using namespace std;

class Environment : public pfate::env::Climate, public pfate::env::LightEnvironment {
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
	C.clim_inst.ppfd = 425;        // umol/m2/s
	C.clim_inst.vpd  = 735;        // Pa
	C.clim_inst.co2  = 368.9;      // ppm
	C.clim_inst.elv = 0;           // m.a.s.l
	C.clim_inst.swp = -0.04;       // MPa

	pfate::env::Clim ca;
	ca.tc = 26.896;         // temperature, deg C
	ca.ppfd = 1866;        // umol/m2/s
	ca.vpd  = 735;        // Pa
	ca.co2  = 368.9;      // ppm
	ca.elv = 0;           // m.a.s.l
	ca.swp = -0.04;       // MPa
	C.set_forcing_acclim(0, ca);

	auto res = P.assimilator.net_production(C, &P.geometry, P.par, P.traits);	

	std::cout << "phydro (whole-plant avg): " << '\n';
	std::cout << "   co = " << res.c_open_avg << '\n'
			  << "   a = " << res.gpp << '\n' 
	          << "   e = " << res.trans << '\n' 
			  << "   vcmax = " << res.vcmax_avg << '\n' 
			  << "   vcmax25 = " << res.vcmax25_avg << '\n'
			  << "   gs = " << res.gs_avg << '\n'
			  << "   dpsi = " << res.dpsi_avg << '\n';

	return 0;
}


