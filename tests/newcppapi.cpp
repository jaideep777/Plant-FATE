#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <netcdf>
using namespace std;
using namespace netCDF;


// print a vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	//os << std::setprecision(12);
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	os << '\n';
	return os;
}


int main(){

	NcFile f("tests/data/gpp.2000-2015.nc", NcFile::read);
	if (!f.isNull() ) cout << "Success\n";
	
	NcVar lonVar = f.getVar("lon");
	if (lonVar.isNull()) cout << "Found lon\n";
	
	std::multimap<std::string, NcDim> dims_map = f.getDims();
	cout << "dims in file: ";
	for (auto p : dims_map) cout << p.first << " "; 
	cout << endl;

	cout << "NcVar lon has size = " << sizeof(NcVar) << " bytes" << endl;
	
	vector <NcDim> dims = lonVar.getDims();
	cout << "Lon var has " << dims.size() << " dims\n";
	
	for (auto dim : dims) cout << dim.getSize() << " ";
	cout << "\n"; 	

	cout << "Lon var has " << lonVar.getDim(0).getSize() << " values\n";

	vector <float> lons(dims[0].getSize());
	lonVar.getVar(lons.data());
	
	string unit;
	lonVar.getAtt("units").getValues(unit);
	cout << "Lon units are: " << unit << endl;
	
	multimap<string,NcVar> vars_map = f.getVars();
	cout << "vars in file: ";
	for (auto p : vars_map) cout << p.first << " "; 
	cout << endl;
	
	map<string,NcGroup> coords_map_temp = f.getCoordVars();
	cout << "coords in file: ";
	for (auto p : coords_map_temp) cout << p.first << " "; 
	cout << endl;

	map<string, NcVar> coords_map;
	for (auto p : coords_map_temp){
		coords_map[p.first] = p.second.getVar(p.first);
	} 
	cout << "coord vars obtained: ";
	for (auto p : coords_map) cout << p.first << " "; 
	cout << endl;

	for (auto p : coords_map){
		vars_map.erase(p.first);
	}	
	cout << "vars left after erasing coords: ";
	for (auto p : vars_map) cout << p.first << " "; 
	cout << endl;

	vector <NcDim> var1_dims = vars_map.begin()->second.getDims();
	vector<std::string> var1_dimnames;
	vector<size_t> var1_dimsizes;
	for (auto dim : var1_dims){
		var1_dimnames.push_back(dim.getName());
		var1_dimsizes.push_back(dim.getSize());
	}

	cout << "First var (" << vars_map.begin()->first << ") has " << var1_dims.size() << " dims: ";
	cout << var1_dimnames;

	cout << "First var (" << vars_map.begin()->first << ") has dimension: ";
	cout << var1_dimsizes;

	// auto it = std::find(var1_dimnames.begin(), var1_dimnames.end(), "time");
	// *it = "time";

	auto itv = std::find(var1_dimnames.begin(), var1_dimnames.end(), "lat");
	*itv = "LAT";

	itv = std::find(var1_dimnames.begin(), var1_dimnames.end(), "lon");
	*itv = "longitude";
	std::cout << "Modified dim names: ";
	for (auto s : var1_dimnames) std::cout << s << " "; 
	std::cout << "\n";

	vector<std::string>   t_names_try = {"time"};
	vector<std::string> lev_names_try = {"lev", "level", "z"};
	vector<std::string> lat_names_try = {"lat", "latitude", "y"};
	vector<std::string> lon_names_try = {"lon", "longitude", "x"};
	map<std::string, std::string> renaming_map;
	for (auto s :   t_names_try) renaming_map[s] = "time";
	for (auto s : lev_names_try) renaming_map[s] = "lev";
	for (auto s : lat_names_try) renaming_map[s] = "lat";
	for (auto s : lon_names_try) renaming_map[s] = "lon";

	// map<std::string, int> coord_ids_map;
	for (auto& name : var1_dimnames){
		std::transform(name.begin(), name.end(), name.begin(),
    	               [](unsigned char c){ return std::tolower(c); });
		if (renaming_map.find(name) != renaming_map.end()){
			name = renaming_map[name];
			// coord_ids_map.insert({name, count})
		}
	}

	std::cout << "Standardized dim names: ";
	for (auto s : var1_dimnames) std::cout << s << " "; 
	std::cout << "\n";

	int unlim_idx = -1;
	for (int i=0; i<var1_dims.size(); ++i) if (var1_dims[i].isUnlimited()) unlim_idx = i; 
	cout << "unlimited / time ID = " << unlim_idx << '\n';

	int lon_idx = std::find(var1_dimnames.begin(), var1_dimnames.end(), "lon") - var1_dimnames.begin();
	int lat_idx = std::find(var1_dimnames.begin(), var1_dimnames.end(), "lat") - var1_dimnames.begin();
	cout << "lat/lon ids = " << lat_idx << " " << lon_idx << "\n";
	if (lon_idx >= var1_dimsizes.size() || lat_idx >= var1_dimsizes.size()) throw std::runtime_error("Lat or Lon not found"); 

	vector<size_t> start(var1_dimnames.size(), 0); // default starts at 0,0,..
	vector<size_t> count = var1_dimsizes;          // default counts as dim sizes (all data will be read at once by default)
	if (unlim_idx != -1) count[unlim_idx] = 1;     // read only 1 time slice
	
	start[lon_idx] = 161;
	start[lat_idx] = 230;
	count[lon_idx] = 1;
	count[lat_idx] = 1;

	vector<double> slice(1);
	vars_map.begin()->second.getVar(start, count, slice.data());
	cout << "slice: " << slice;

	// cout << "Lon var has " << lonVar.getDim(0).getSize() << " values\n";


//	for (auto it : coords_map) cout << it.first << " ";


	
	// Get the coordinate variables
	NcFile* dFile = &f;
	NcVar latVar, levVar, tVar;
	int ncoords = 0;
	for (auto p : coords_map){
		auto it = &p;

		string var = it->first;
		
		if (it->first == "lat" || it->first == "latitude" || it->first == "LAT" ){
			latVar = it->second;
			++ncoords;
		}
//		else cout << "lat not found!\n";

		if (it->first == "lon" || it->first == "longitude" || it->first == "LON" ){
			lonVar = it->second;
			++ncoords;
		}
//		else cout << "lon not found!\n";


		if (it->first == "lev" || it->first == "levels" || it->first == "LEV" ){
			levVar = it->second;
			++ncoords;
		}
//		else cout << "lev not found!\n";


		if (it->first == "time" || it->first == "TIME" || it->first == "LEV" ){
			tVar = it->second;
			++ncoords;
		}
//		else cout << "time not found!\n";

	}
	
	if (latVar.isNull()) cout << "Lat not found\n";
	if (lonVar.isNull()) cout << "Lon not found\n";
	if (levVar.isNull()) cout << "Lev not found\n";
	if (tVar.isNull())   cout << "Time not found\n";
	
	
	NcVar firstVar;
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it ){
		if (it->second.getDimCount() == ncoords) {firstVar = it->second; break;}
	}
	cout << firstVar.getName() << " (" << firstVar.getId() << ")";
	
	string a;
	NcVarAtt A = firstVar.getAtt("units");
	if (!A.isNull()) A.getValues(a);
	else a = "NA";
	cout << "var unit = " << a << endl;
	
	// missing value
	float b;
	try{
		A = firstVar.getAtt("missing_value");
		if (!A.isNull()) A.getValues(&b);
	}
	catch(...){
		cout << "missing_value not found\n";
	}

	try{
		A = firstVar.getAtt("_FillValue");
		if (!A.isNull()) A.getValues(&b);
	}
	catch(...){
		cout << "_FillValue not found\n";
	}
	cout << "missing Value = " << b << endl;

	
	return 0;
} 
