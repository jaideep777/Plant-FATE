#ifndef PLANT_FATE_ENV_CLIMATE_H_
#define PLANT_FATE_ENV_CLIMATE_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>


namespace env{

class Clim{
	public:
	double tc = 25.5;         // temperature, deg C
	double ppfd_max = 1000;   // umol/m2/s
	double ppfd = 377;        // umol/m2/s
	double vpd  = 540;        // Pa
	double co2  = 368.9;      // ppm
	double elv = 0;           // m.a.s.l
	double swp = -0.04;       // MPa

};


class Climate{
	
	private:
	double t_prev = 0;  // time at current data values (years since 2000-01-01)
	double t_next = 0;  // next time in file for which data is available
	Clim clim_prev;
	Clim clim_next;
	double t0 = 2000.0;
	double t_base = 2000.0;
	double delta = 0;

	public:
	double t_now;
	Clim clim;

	std::string metFile = "";
	std::string co2File = "";
	bool interpolate = false;

	private:
	std::ifstream fin_met;
	std::ifstream fin_co2;

	public:
	
	int init();
	
	Clim interp(Clim &clim_prev, Clim &clim_next);

	void updateClimate(double t);

	int readNextLine_met();

	template<class T> 
	T as(std::string s);

	void print(double t);

};


template<class T> 
T Climate::as(std::string s){
	std::stringstream sin(s);
	T data;
	sin >> data;
	return data;
}


} // namespace env


#endif
