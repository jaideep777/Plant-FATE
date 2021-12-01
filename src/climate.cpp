
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

#include "climate.h"


namespace env{

		
int Climate::init(){
	fin_met.open(metFile.c_str());
	fin_co2.open(co2File.c_str());
	
	if (!fin_met){
		std::cerr << "Error: Could not open file " << metFile << "\n";
		return 1;
	}
	if (!fin_co2){
		std::cerr << "Error: Could not open file " << co2File << "\n";
		return 1;
	}
	
	// skip header
	std::string line;
	getline(fin_met, line);
	getline(fin_co2, line);

	// read first entry into t_prev
	readNextLine_met();
	t_now = t_prev = t_next;
	clim = clim_prev = clim_next;
	readNextLine_met();
	
	return 0;

}


Clim Climate::interp(Clim &clim_prev, Clim &clim_next){
	return clim_prev;
}

void Climate::updateClimate(double t){
	if (t == t_now) return;
	
	if (t > t_prev && t < t_next){
		clim = interp(clim_prev, clim_next);
	}

	while (t >= t_next){
		std::cout << "update - " << int(t) << "/" << (t-int(t))*12+1 << " --> " << int(t_next) << "/" << (t_next-int(t_next))*12+1 <<  "\n"; std::cout.flush();
		clim_prev = clim_next;
		t_prev = t_next;
		readNextLine_met();
	}

	clim = interp(clim_prev, clim_next);		
}


int Climate::readNextLine_met(){

	//std::vector<std::string>   result;
	std::string                line, cell;
	
	if (fin_met.peek() == EOF){
		std::cout << "RESET FILE\n"; std::cout.flush();
		fin_met.clear();
		fin_met.seekg(0);
		std::getline(fin_met, line); // skip header
		t_base = t_next+1/12.0;
	}
	
	// READ line
	std::getline(fin_met, line);

	std::stringstream          lineStream(line);
	
	// PARSE line
	std::getline(lineStream, cell, ',');
	int year = as<int>(cell);
	
	std::getline(lineStream, cell, ',');
	int month = as<int>(cell);

	std::getline(lineStream, cell, ',');
	clim_next.tc = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim_next.vpd = as<double>(cell)*100;  // convert hPa to Pa
	
	std::getline(lineStream, cell, ',');
	clim_next.ppfd = as<double>(cell);

	std::getline(lineStream, cell, ',');
	clim_next.ppfd_max = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim_next.swp = as<double>(cell);

	t_next = t_base + (year-t0) + (month-1)/12.0;
		
	return 0;
}


void Climate::print(double t){
	int year = int(t);
	double month = (t-int(t))*12+1;
	year += int(month/12.0); month = fmod(month, 12.0);
	int yearp = int(t_prev);
	double monthp = (t_prev-int(t_prev))*12+1;
	int yearn = int(t_next);
	double monthn = (t_next-int(t_next))*12+1;
	std::cout << "Climate at t = " << year << "/" << month << "\n";
	std::cout << "prev: " << yearp << "/" << monthp << " | " << clim_prev.vpd << " " << clim_prev.ppfd << "\n"; 
	std::cout << "now : " << year  << "/" << month  << " | " << clim.vpd      << " " << clim.ppfd      << "\n"; 
	std::cout << "next: " << yearn << "/" << monthn << " | " << clim_next.vpd << " " << clim_next.ppfd << "\n"; 
}

} // namespace env


