
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"


namespace env{


int Climate::init(){
	fin_met.open(metFile.c_str());
	fin_co2.open(co2File.c_str());
	
	if (!fin_met){
		throw std::runtime_error("Could not open file " + metFile);
	}
	if (!fin_co2){
		throw std::runtime_error("Could not open file " + co2File);
	}
	
	// skip header
	std::string line;
	getline(fin_met, line);
	getline(fin_co2, line);

	while (fin_met.peek() != EOF){
		Clim clim1;
		double t1;
		readNextLine_met(clim1, t1);
		t_met.push_back(t1);
		v_met.push_back(clim1);
	}
	
	delta = *t_met.rbegin() - *t_met.begin() + t_met[1]-t_met[0]+ 1e-6;
	
	// read CO2 file
	while (fin_co2.peek() != EOF){
		std::getline(fin_co2, line);
	
		std::stringstream lineStream(line);

		std::string cell;
		
		std::getline(lineStream, cell, ',');
		int year = as<int>(cell);

		std::getline(lineStream, cell, ',');
		double co2 = as<double>(cell);
		
		t_co2.push_back(year);
		v_co2.push_back(co2);
	}
	
	return 0;

}


Clim Climate::interp(Clim &clim_prev, Clim &clim_next){
	return clim_prev;
}

int Climate::id(double t){
	int id = (t - *t_met.begin())/delta*t_met.size();
	id = id % t_met.size();
	return id;
}


void Climate::updateClimate(double t){
		
	if (update_met){
		double tadj = t;  // adjusted t to lie between limits of observed data
		while(tadj < *t_met.begin()) tadj += delta;
		int idx_now = id(tadj);
		int idx_next = (idx_now+1) % t_met.size();
		clim = interp(v_met[idx_now], v_met[idx_next]);
	}

	if (update_co2){
		int year = int(t);
		//std::cout << "CO2: " << t << " " << *t_co2.begin() << " " << *t_co2.rbegin() << "\n";
		if (year >=  *t_co2.begin() && year <= *t_co2.rbegin()){
			clim.co2 = v_co2[year - t_co2[0]];
			//std::cout << "Setting co2 @ t = " << year << " from (" << year - t_co2[0] << " / " << v_co2[year - t_co2[0]] << ")\n";
		}
		else if (year > *t_co2.rbegin()) {
			//std::cout << "Setting co2 @ end \n";
			clim.co2 = *v_co2.rbegin();
		}
	}

}


int Climate::readNextLine_met(Clim &clim, double &t){

	std::string                line, cell;
	
	// READ line
	std::getline(fin_met, line);

	std::stringstream          lineStream(line);
	
	// PARSE line
	std::getline(lineStream, cell, ',');
	int year = as<int>(cell);
	
	std::getline(lineStream, cell, ',');
	int month = as<int>(cell);

	std::getline(lineStream, cell, ',');
	clim.tc = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim.vpd = as<double>(cell)*100;  // convert hPa to Pa
	
	std::getline(lineStream, cell, ',');
	clim.ppfd = as<double>(cell);

	std::getline(lineStream, cell, ',');
	clim.ppfd_max = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim.swp = -as<double>(cell);  // convert to negative (in file it is absolute value)

	t = year + (month-1)/12.0;
		
	return 0;
}


void Climate::print(double t){
	int year = int(t);
	double month = (t-int(t))*12;
	std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
	std::cout << " | " << clim.tc << " " << clim.ppfd      << " " << clim.vpd << " " << clim.co2 << "\n"; 
}


void Climate::print_all(){
	for (int i = 0; i<t_met.size(); ++i){
		double t = t_met[i];
		int year = int(t);
		double month = (t-int(t))*12;
		std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
		std::cout << " | " << v_met[i].tc << " " << v_met[i].ppfd      << " " << v_met[i].vpd << "\n"; 
	}
}


} // namespace env


