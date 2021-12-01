#ifndef PLANT_FATE_IO_INITIALIZER_H_
#define PLANT_FATE_IO_INITIALIZER_H_

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
/** \ingroup utils */

/**
	\brief A simple initializer that reads a parameter file and stores the values in a named std::map.
	
	The parameter file must have three sections - STRINGS, SCALARS, ARRAYS. Sections start with '>'. 
	Each section has name-value pairs separated by whitespace. Arrays have a name followed by values, ending in '-1'.
	Comments are allowed. Comments start with "# " (note the space) and can come either on a new line or on the same line after
	the name-value(s) pair. 
	
	Here is an example parameter file:
	
	~~~{.ini}
	> STRINGS
	sim_name      mySimution
	output_file   ~/output/test.txt
	
	> SCALARS
	graphics      1           # Do we want graphics to be on? 
	timesteps     1000        # For how many timesteps do we run the simulation?
	dt            0.1
	# This is also a comment
	
	> ARRAYS
	parameter1    1 2 3 4 5 6 -1
	~~~
*/
namespace io{

class Initializer{
	private:
	std::string init_fname;
	std::ifstream fin;
	
	std::map <std::string, std::string> strings;
	std::map <std::string, double>  scalars;
	std::map <std::string, std::vector <double>> arrays;
	
	public:
	
	inline Initializer(){}

	inline Initializer(std::string fname){
		init_fname = fname;
	}

	inline void setInitFile(std::string fname){
		init_fname = fname;
	}

	inline int readFile(){
		// Reset maps here
		strings.clear(); 
		scalars.clear();
		arrays.clear();

		fin.open(init_fname.c_str());
		if (!fin) {
			std::cerr << "FATAL ERROR: Cannot open initializer file " << init_fname << '\n';
			return 1;
		}
		
		std::string attr_begin = ">";
		std::string init_format = "    > STRINGS\n    ... \n\n    > SCALARS\n    ...\n\n    > ARRAYS \n    ...\n";

		std::string s, v;
		double f;
		std::vector <double> vf;
		
		while (fin >> s && s != attr_begin);	// read until 1st > is reached

		fin >> s; 
		if (s != "STRINGS") {
			std::cout << "FATAL ERROR: STRINGS section missing in Initializer file. Format must be:\n" << init_format << '\n'; 
			return 1;
		}
		while (fin >> s && s != attr_begin){
			if (s == "") continue;	// skip empty lines
			if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #-followed lines (comments)
			fin >> v;
			strings[s] = v;
		}

		fin >> s;
		if (s != "SCALARS") {
			std::cout << "FATAL ERROR: SCALARS section missing in Initializer file. Format must be:\n" << init_format << '\n'; 
			return 1;
		}
		while (fin >> s && s != attr_begin){
			if (s == "") continue;	// skip empty lines
			if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #-followed lines (comments)
			fin >> f;
			scalars[s] = f;
		}

		fin >> s;
		if (s != "ARRAYS") {
			std::cout << "FATAL ERROR: ARRAYS section missing in Initializer file. Format must be:\n" << init_format << '\n'; 
			return 1;
		}
		while (fin >> s && s != attr_begin){
			if (s == "") continue;	// skip empty lines
			if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #-followed lines (comments)
			while (fin >> f && f != -1) vf.push_back(f);	// TODO: Replace -1 with a symbol
			arrays[s] = vf;
			vf.resize(0);
		}

		fin.close();
	}

	template<class T>
	T get(std::string s){
		std::map <std::string, std::string>::iterator it = strings.find(s);
		if (it != strings.end()){
			std::string val_s = it->second;
			std::stringstream sin(val_s);
			T val;
			sin >> val;
			return val;	
		}
		else {
			std::cout << "FATAL ERROR: Could not find required variable " << s << " in initializer file.\n";
			return T();
		}
		
	}
	
	
	inline std::string getString(std::string s){
		std::map <std::string, std::string>::iterator it = strings.find(s);
		if (it != strings.end()) return it->second;
		else {
			std::cout << "FATAL ERROR: Could not find required variable " << s << " in initializer file.\n";
			return "";
		}
	}

	inline double getScalar(std::string s){
		std::map <std::string, double>::iterator it = scalars.find(s);
		if (it != scalars.end()) return it->second;
		else {
			std::cout << "FATAL ERROR: Could not find required variable " << s << " in initializer file.\n";
			return 1;
		}
	}

	inline std::vector <double> getArray(std::string s, int size = -1){
		std::map <std::string, std::vector<double> >::iterator it = arrays.find(s);
		if (it == arrays.end()) {	// array not found
			std::cout << "FATAL ERROR: Could not find required array " << s << " in initializer file.\n";
			return {};
		}
		if (it->second.size() == 0){
			std::cout << "FATAL ERROR: Required array " << s << " is empty!\n"; 
			return {};
		}
		if (size == -1) return it->second;
		else if (size != it->second.size()) {
			std::cout << "FATAL ERROR: Incorrect size of array " << s << ". Required " << size << ", found " << it->second.size() << '\n'; 
			return {};
		}
		else return it->second;
	}

	inline void print(){
		std::cout << "-------:\n";
		std::cout << "STRINGS:\n";
		std::cout << "-------:\n";
		for (std::map<std::string, std::string>::iterator it = strings.begin(); it != strings.end(); ++it){
			std::cout << it->first << ": " << it->second << '\n';
		}
		std::cout << "-------:\n";
		std::cout << "SCALARS:\n";
		std::cout << "-------:\n";
		for (std::map<std::string, double>::iterator it = scalars.begin(); it != scalars.end(); ++it){
			std::cout << it->first << ": " << it->second << '\n';
		}
		std::cout << "-------:\n";
		std::cout << "ARRAYS:\n";
		std::cout << "------:\n";
		for (std::map<std::string, std::vector<double> >::iterator it = arrays.begin(); it != arrays.end(); ++it){
			std::cout << it->first << ": ";
			for (int i=0; i< it->second.size(); ++i) std::cout << it->second[i] << " ";
			std::cout << "\n";
		}
		std::cout << "\n";
	}
	
	
};

} // namespace io

#endif

