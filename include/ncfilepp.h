#ifndef FLARE_FLARE_NCFILEPP_H
#define FLARE_FLARE_NCFILEPP_H

#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <netcdf>
#include <chrono>
#include <cmath>
#include <limits>

#include "utils.h"

namespace flare{

class CoordMeta{
	public:
	std::string name;
	std::string unit;
	std::vector<double> values;	
};


class VarMeta{
	public:
	std::string name;
	std::string unit;

	std::vector <std::string> dimnames;
	std::vector <size_t> dimsizes;
	 
	int unlim_idx = -1;
};


class NcFilePP : public netCDF::NcFile {
	public:
	std::map <std::string, netCDF::NcVar> vars_map;     // variable name --> NcVar map (for data variables)
	std::map <std::string, netCDF::NcVar> coords_map;   // variable name --> NcVar map (for coordinate variables)

	std::map <std::string, CoordMeta> coords; 
	std::map <std::string, VarMeta> vars; 

	// data for standardizing dimension names (kept public to be editable)
	std::vector<std::string>   t_names_try = {"time", "t"};
	std::vector<std::string> lev_names_try = {"lev", "level", "z"};
	std::vector<std::string> lat_names_try = {"lat", "latitude", "y"};
	std::vector<std::string> lon_names_try = {"lon", "long", "longitude", "x"};

	std::map<std::string, std::string> renaming_map; // map used for renaming dimensions to standard names

	std::string standarize_name(std::string s){
		// convert name to lowercase (all keys in renaming map are in lowercase)
		std::transform(s.begin(), s.end(), s.begin(),
		               [](unsigned char c){ return std::tolower(c); });

		// check if it exists in renaming map
		auto it = renaming_map.find(s);

		if (it != renaming_map.end()) return it->second; // if entry is found, return the standardised name
		else return s;                                   // else return original string (converted to lowercase)
	}

	inline void readMeta(){
		// Create renaming map for dimension names
		for (auto s :   t_names_try) renaming_map[s] = "time";
		for (auto s : lev_names_try) renaming_map[s] = "lev";
		for (auto s : lat_names_try) renaming_map[s] = "lat";
		for (auto s : lon_names_try) renaming_map[s] = "lon";

		// get all variables in the file in a [name --> variable] map
		std::multimap<std::string, netCDF::NcVar> vars_map_temp = this->getVars();
		// then create vars_map by standardizing names
		for (auto& p : vars_map_temp){
			vars_map[standarize_name(p.first)] = p.second;
		}
		
		// Get all coordinates. This returns NcGroups
		std::map<std::string, netCDF::NcGroup> coords_map_temp = this->getCoordVars();
		// Extract variables from the obtained coordinate names
		for (auto& p : coords_map_temp){
			coords_map[standarize_name(p.first)] = p.second.getVar(p.first);
		} 

		// remove coords from list of all variables, so only data variables are left
		for (auto p : coords_map){
			vars_map.erase(p.first);
		}

		// read values and units for coordinate vars
		for (auto p : coords_map){
			CoordMeta c;
			c.name = standarize_name(p.first);

			std::string s;
			try{
				p.second.getAtt("units").getValues(s);
				c.unit = s;
			}
			catch(...){
				std::cout << "Warning: variable " << p.first << " has no unit.\n";
				c.unit = "";
			}

			c.values.resize(p.second.getDim(0).getSize());
			p.second.getVar(c.values.data());

			coords[standarize_name(p.first)] = c;
		}

		// read metadata for data variables
		for (auto p : vars_map){
			VarMeta vm;
			netCDF::NcVar& ncvar = p.second;
			vm.name = ncvar.getName();

			std::vector <netCDF::NcDim> ncdims = ncvar.getDims();

			// get index of unlimited dimension
			for (int i=0; i<ncdims.size(); ++i){
				if (ncdims[i].isUnlimited()) vm.unlim_idx = i;
			}

			// get names and sizes of dimensions for this variable
			for (auto d : ncdims){
				vm.dimnames.push_back(standarize_name(d.getName()));
				vm.dimsizes.push_back(d.getSize());
			}

			// Get basic variable attributes:
			// unit
			try{ ncvar.getAtt("units").getValues(vm.unit); }
			catch(netCDF::exceptions::NcException &e){ std::cout << "Warning: Variable does not have a unit\n";}

			vars[vm.name] = vm;
		}

	}

	inline void printMeta(){
		std::cout << "\nNC file >" << "\n";
		std::cout << "   coords:\n";
		for (auto& p : coords) std::cout << "      " << p.first << " (" << p.second.unit << ")\n";
		std::cout << "   vars (excluding coords):\n";
		for (auto p : vars){
			auto& v = p.second;
			std::cout << "      " << v.name << " (" << v.unit << ")\n";
			std::cout << "         dims: " << v.dimnames << "\n";
			std::cout << "         dim sizes: " << v.dimsizes << "\n";
			std::cout << "         unlimited dim: " << ((v.unlim_idx >= 0)? v.dimnames[v.unlim_idx] : "NA") << "\n";
		}
		std::cout << "   coord values:\n";
		for (auto p : coords){
			std::cout << "      " << p.first << ": ";
			if (p.second.values.size() <= 6) std::cout << p.second.values << '\n';
			else{
				auto& v = p.second.values;
				size_t n = p.second.values.size();
				std::cout << n << " | " << v[0] << " " << v[1] << " " << v[2] << " ... " 
				          << v[n-3] << " " << v[n-2] << " " << v[n-1] << "\n";
			}
		}
		std::cout << "~~\n";
	}
};


} // namespace flare

#endif
