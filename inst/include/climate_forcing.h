#ifndef PLANT_FATE_ENV_CLIMATE_FORCING_H_
#define PLANT_FATE_ENV_CLIMATE_FORCING_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "climate.h"


namespace env{

class ClimateForcing{
	private:
	double t_prev = 0;  // time at current data values (years since 2000-01-01)
	double t_next = 0;  // next time in file for which data is available
	Clim clim_prev;
	Clim clim_next;
   	double t0 = 2000.0;
	double tf = 2001 + 11.0/12;
	double t_base = 2000.0;
	double delta = 1.0/12;
    double deltaT = 1;

	std::vector<double>	t_co2;
	std::vector<double> v_co2;

	//std::vector<double> t_met;
	//std::vector<Clim>   v_met;
	// Adding temp vector for soil water potential
	//std::vector<double> v_met_swp;

	public:
	double t_now = 0;
	Clim clim;
	std::vector<double> t_met;
	std::vector<Clim>   v_met;

	std::string metFile = "";
	std::string co2File = "";
	bool interpolate = false;
	
	bool update_met = false;
	bool update_co2 = false;

	private:
	// std::ifstream fin_met; // Cannot use streams as members because when exported to R, R needs copy constructor
	// std::ifstream fin_co2;

	public:
	
	int init();
	
	void set(double _tc, double _ppfd_max, double _ppfd, double _vpd, double _co2, double _elv, double _swp);

	void update_tc(double t, double _tc);
	void update_ppfd_max(double t, double _ppfd_max);
	void update_ppfd(double t, double _ppfd);
	void update_vpd(double t, double _vpd);
	// void update_co2(double t, double _co2);
	void update_elv(double t, double _elv);
	void update_swp(double t, double _swp);

	Clim interp(Clim &clim_prev, Clim &clim_next);

	int id(double t); // since I'm not sure what this is doing I'm adding the idx function here
    int idx_of(double t);
	void updateClimate(double t);

	int readNextLine_met(Clim &clim, double &t, std::ifstream& fin_met);
	int readNextLine_co2(double &co2, double &t, std::ifstream& fin_co2);
	
	template<class T> 
	T as(std::string s);

	void print();
	void print_line(double t);
	void print_all();

	private: 
	double findMinDeltaT(std::vector <double> t_vector);

};


template<class T> 
T ClimateForcing::as(std::string s){
	std::stringstream sin(s);
	T data;
	sin >> data;
	return data;
}


} // namespace env


#endif
