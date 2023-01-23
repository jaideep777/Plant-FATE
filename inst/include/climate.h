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


class Climate{
	private:
	double t_prev = 0;  // time at current data values (years since 2000-01-01)
	double t_next = 0;  // next time in file for which data is available
	Clim clim_prev;
	Clim clim_next;
   	double t0 = 2000.0;
	double tf = 2001 + 11.0/12;
	double t_base = 2000.0;
	double delta = 1.0/12;

	std::vector<int>    t_co2;
	std::vector<double> v_co2;

	//std::vector<double> t_met;
	//std::vector<Clim>   v_met;
	// Adding temp vector for soil water potential
	//std::vector<double> v_met_swp;

	public:
	double t_now;
	Clim clim;
	std::vector<double> t_met;
	std::vector<Clim>   v_met;
	// Adding temp vector for soil water potential
	std::vector<double> v_met_swp;
	int counter_var = 1;

	std::string metFile = "";
	std::string co2File = "";
	bool interpolate = false;
	
	bool update_met = true;
	bool update_co2 = true;

	private:
	// std::ifstream fin_met; // Cannot use streams as members because when exported to R, R needs copy constructor
	// std::ifstream fin_co2;

	public:
	
	int init();
	
	Clim interp(Clim &clim_prev, Clim &clim_next);

	int id(double t);
	void updateClimate(double t);

	int readNextLine_met(Clim &clim, double &t, std::ifstream& fin_met);
	
	int binarySearch(double k);
	double inst_swp(double year);
	
	template<class T> 
	T as(std::string s);

	void print(double t);
	void print_all();

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
