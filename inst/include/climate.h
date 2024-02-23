#ifndef PLANT_FATE_ENV_CLIMATE_H_
#define PLANT_FATE_ENV_CLIMATE_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <phydro.h>  // for calc_patm()

namespace env{

class Clim{
	public:
	double tc = 25.5;         // Temperature [C]
	double ppfd_max = 2000;   // PAR (daytime max) [umol m-2 s-1]
	double ppfd = 500;        // PAR (daily 24-hr mean) [umol m-2 s-1]
	double rn = 250;          // Net radiation at surface [W m-2]
	double vpd  = 540;        // Vapour pressure deficit [Pa]
	double co2  = 368.9;      // Atmospheric CO2 [ppm]
	double elv = 0;           // Site elevation [m.a.s.l]
	double swp = -0.04;       // Soil water potential [MPa]
	double pa;                // Surface pressure [Pa]

	Clim(){
		pa = phydro::calc_patm(elv);
	}

	void set_elevation(double _elv){
		elv = _elv;
		pa = phydro::calc_patm(elv);
	}
};


// TODO: for now, this only has one clim object, but this class 
// should store instantaneous and acclimation forcing, perform
// running averages, etc.
class Climate{
	public:
	Clim clim_inst;
	Clim clim_acclim;

	void set_elevation(double _elv);
	void set_co2(double _co2);

	virtual void print(double t);
	void print_line(double t);
};

} // namespace env


#endif
