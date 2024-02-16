#ifndef PLANT_FATE_ENV_CLIMATE_H_
#define PLANT_FATE_ENV_CLIMATE_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>


namespace env{

class Clim{
	public:
	double tc = 25.5;         // temperature, deg C
	double ppfd_max = 2000;   // umol/m2/s
	double ppfd = 500;        // umol/m2/s
	double vpd  = 540;        // Pa
	double co2  = 368.9;      // ppm
	double elv = 0;           // m.a.s.l
	double swp = -0.04;       // MPa

};


// TODO: for now, this only has one clim object, but this class 
// should store instantaneous and acclimation forcing, perform
// running averages, etc.
class Climate{
	public:
	Clim clim;

	virtual void print(double t);
	void print_line(double t);
};

} // namespace env


#endif
