#ifndef PLANT_FATE_IO_TRAITS_READER_H_
#define PLANT_FATE_IO_TRAITS_READER_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>

#include "traits_params.h"

class TraitsReader{
	public:
	
	std::vector<plant::PlantTraits> species;
	
	template<class T> 
	T as(std::string s){
		std::stringstream sin(s);
		T data;
		sin >> data;
		return data;
	}
	
	inline int readFromFile(std::string fname){

		std::ifstream fin(fname.c_str());
		if (!fin){
			throw std::runtime_error("Could not open file " + fname + "\n");
		}

		//std::vector<std::string>   result;
		std::string line, cell;
		getline(fin, line); // ignore header
		
		while (fin.peek() != EOF){
			// READ line
			std::getline(fin, line);

			std::stringstream          lineStream(line);
			
			plant::PlantTraits traits;

			// PARSE line
			std::getline(lineStream, cell, ','); // family
			
			std::getline(lineStream, cell, ','); // species name
			traits.species_name = cell;	

			std::getline(lineStream, cell, ','); // num ind
			std::getline(lineStream, cell, ','); // basal area
			
			std::getline(lineStream, cell, ','); // wd
			if (cell != "" && cell != "NA")	traits.wood_density = as<double>(cell)*1000; // g/cc to kg/m3
			else traits.wood_density = 686.638;

			std::getline(lineStream, cell, ','); // hmat
			if (cell != "" && cell != "NA")	traits.hmat = as<double>(cell);
			else traits.hmat = 23.99;

			std::getline(lineStream, cell, ','); // lma
			if (cell != "" && cell != "NA")	traits.lma = as<double>(cell)*1e-3;  // convert g/m2 to kg/m2
			else traits.lma = 0.119378;

			std::getline(lineStream, cell, ','); // p50
			if (cell != "" && cell != "NA")	traits.p50_xylem = as<double>(cell); 
			else traits.p50_xylem = -2.29;

			// ignore further data (for now)			

			species.push_back(traits);
		}
		
		return 0;
	}

	inline void print(){
		for (auto& s : species){
			std::cout << s.species_name << ":  " << s.lma << "\t" << s.wood_density << "\t" << s.hmat << "\t" << s.p50_xylem << "\n";
		}

	}

};



#endif


