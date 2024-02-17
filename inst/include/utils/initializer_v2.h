#ifndef UTILS_IO_INITIALIZER_H_
#define UTILS_IO_INITIALIZER_H_

#include <iostream>
#include <regex>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

/**
	\brief A simple initializer that reads parameters from an ini file.
	
	reference: https://codereview.stackexchange.com/questions/127819/ini-file-parser-in-c
	
	The parameter file must follow the formatting requirements of .ini files. 
	Sections are enclosed with [] and cannot have spaces or comments on the 
	same line. Each section has name-value pairs separated by =. Arrays have 
	values seperated by whitespace. Comments start with a ";" or "#" 
	and can occur anywhere in the file except on section lines. 
	All content on a line folling the 
	comment character will be ignored. 
	
	Here is an example parameter file:
	
	~~~{.ini}
	; a comment
	; this unnamed section is treated as the global section
	sim_name    =   mySimution    ;another comment
	output_file =  ~/output/test.txt  # yet another comment
	
	# this is also a comment
	; comments cannot be inserted in combinaiton with section headers. 
	; Also, no spaces after ] below
	[section 2]
	graphics    =  1           # Do we want graphics to be on? 
	timesteps   =  1000        ; For how many timesteps do we run the simulation?
	dt          =  0.1
	
	[arrays]
	array1      = 1 2 3 4 5 6  # these numbers can be retreived in a vector
	~~~

	Use:
	----
	io::Initializer I;
	I.parse(file);
	
	I.get<double>("sim_name");  
	I.get<double>("section 2", "dt");  
	I.get_vector<double>("array", "array1");

*/
namespace io{

class Initializer{
	using section = std::unordered_map<std::string, std::string>;

	private:
	std::unordered_map<std::string, section> sections;
	std::ifstream fin;

	private:
	inline const section& get_section(const std::string& sectionname) const {
		auto found = sections.find(sectionname);
		if (found != sections.end()) return found->second;
		else throw std::runtime_error("Initializer: Cannot find required section ["+sectionname+"]");
	}
	
	inline std::string get_value(const std::string& sectionname, const std::string& keyname) const {
		section sect = get_section(sectionname);
		auto it = sect.find(keyname);
		if (it != sect.end()) return it->second;
		else throw std::runtime_error("Initializer: Could not find required variable [" + keyname + "] in section [" + sectionname + "]");
	}

	public:
	inline void parse(std::istream& in, bool add = false, bool verbose = false){
		if (!add) sections.clear();
		
//		static const std::regex comment_regex{R"x(\s*[;#])x"};
		static const std::regex section_regex{R"(\s*\[([^\]]+)\])"};
		static const std::regex value_regex{R"(\s*(\S[^ \t=]*)\s*=\s*((\s*\S+)+)\s*$)"};
		static const std::regex comment_regex{"([^;#]*)([;#])"};
		std::string current_section = "global";
		std::smatch pieces;
		std::string line;
		while (std::getline(in, line)){
			// trim text following comment characters
			std::regex_search(line, pieces, comment_regex);
			if (pieces.size() == 3){
				if (verbose) std::cout << "Trimming comment line [" << line << "] to [" << pieces[1].str() << "]\n";
				line = pieces[1].str();
			} 
			
			// parse line
			if (line.empty()) {
				// skip comment lines and blank lines					
			}
			else if (std::regex_match(line, pieces, section_regex)) {
				if (pieces.size() == 2) { // exactly one match
					current_section = pieces[1].str();
					if (verbose) std::cout << "--- Section = " << current_section << " ---\n";
				}
			}
			else if (std::regex_match(line, pieces, value_regex)) {
				if (pieces.size() == 4) { // exactly enough matches
					sections[current_section][pieces[1].str()] = pieces[2].str();
					if (verbose) std::cout << pieces[1].str() << " = " << pieces[2].str() << "\n";
				}
			}
			else {
				if (verbose) std::cout << "skipping line [" << line << "]\n";
				//throw std::runtime_error("Cannot parse line "+line);
			}
		}	
	}
		
	inline void parse(std::string filename, bool add = false, bool verbose = false) {
		fin.open(filename.c_str());
		if (!fin) throw std::invalid_argument("Initializer: Could not open file: "+filename);
		if (verbose) std::cout << "Parsing file: " << filename << "\n";
		parse(fin, add, verbose);
	}

	template<class T>
	T get(const std::string& sectionname, const std::string& keyname) const {
		std::string result = get_value(sectionname, keyname);
		std::stringstream sin(result);
		T val;
		sin >> val;
		return val;	
	}

	template<class T>
	T get(const std::string& keyname) const {
		return get<T>("global", keyname);
	}

	template<class T>
	std::vector<T> get_vector(const std::string& sectionname, const std::string& keyname) const {
		std::string result = get_value(sectionname, keyname);
		std::stringstream sin(result);
		T val;
		std::vector<T> vec;
		while(sin >> val){
			vec.push_back(val);
		}
		return vec;	
	}

	template<class T>
	std::vector<T> get_vector(const std::string& keyname) const {
		return get_vector<T>("global", keyname);
	}
	
	inline void print() const {
		std::cout << "------------------\n";
		std::cout << "> Initializer:\n";
		for (const auto& sec : sections){
			std::cout << "  [" << sec.first << "]\n";
			for (const auto& x : sec.second){
				std::cout << "    " << x.first << " = {" << x.second << "}\n";
			}
		}
		std::cout << "------------------\n";
	}

};


} // namespace io

#endif

