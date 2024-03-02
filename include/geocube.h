#ifndef FLARE_FLARE_GEOCUBE_H
#define FLARE_FLARE_GEOCUBE_H

#include <tensor.h>
#include "ncstream.h"
#include "time_math.h"

namespace flare{

template <class T>
class GeoCube : public Tensor<T> {
	public:
	VarMeta meta;

	std::vector<std::vector<double>> coords, coords_trimmed;

	int lon_idx, lat_idx, t_idx;

	T scale_factor = 1.0, add_offset = 0.0; 

	private:
	netCDF::NcVar ncvar;

	std::vector <size_t> starts, counts;
	std::vector <ptrdiff_t> strides;

	bool streaming = false;
	int current_file_id = -999;

	public:

	void readMeta(NcFilePP &in_file, std::string varname = ""){
		// get the named variable, or the first variable in the file
		if (varname != ""){
			ncvar = in_file.vars_map.at(varname);
		}
		else{
			varname = in_file.vars_map.begin()->first;
			ncvar = in_file.vars_map.begin()->second;
		}

		meta = in_file.vars.at(varname);
		for (auto& d : meta.dimnames){
			coords.push_back(in_file.coords.at(d).values);
		}
		coords_trimmed = coords;

		// Get basic variable attributes
		// missing value
		try{ ncvar.getAtt("missing_value").getValues(&this->missing_value);}
		catch(netCDF::exceptions::NcException &e1){
			try{ ncvar.getAtt("_FillValue").getValues(&this->missing_value);}
			catch(netCDF::exceptions::NcException &e2){ std::cout << "Missing/Fill value not found. Setting to NaN\n";}
		}

		// unit
		try{ ncvar.getAtt("units").getValues(meta.unit); }
		catch(netCDF::exceptions::NcException &e){ std::cout << "Warning: Variable does not have a unit\n";}

		// scale factor
		try{ ncvar.getAtt("scale_factor").getValues(&scale_factor); }
		catch(netCDF::exceptions::NcException &e){ scale_factor = 1.f;}
		
		// offset
		try{ ncvar.getAtt("add_offset").getValues(&add_offset);}
		catch(netCDF::exceptions::NcException &e){ add_offset = 0.f;}

		// ~~ Get the dimension indices. i.e., which index in ncdim std::vector is lat, lon, etc
		lon_idx = std::find(meta.dimnames.begin(), meta.dimnames.end(), "lon")  - meta.dimnames.begin();
		lat_idx = std::find(meta.dimnames.begin(), meta.dimnames.end(), "lat")  - meta.dimnames.begin();
		t_idx   = std::find(meta.dimnames.begin(), meta.dimnames.end(), "time") - meta.dimnames.begin();

		if (lon_idx >= meta.dimsizes.size() || lat_idx >= meta.dimsizes.size()) throw std::runtime_error("Lat or Lon not found"); 
		if (t_idx >= meta.dimnames.size()) t_idx = -1;

		// by default, set to read entire geographic extent at 1 frame in the unlimited dimension 
		starts.clear(); starts.resize(meta.dimnames.size(), 0);   // set all starts to 0
		counts = meta.dimsizes;                                   // set all counts to dimension size (all elements)
		if (meta.unlim_idx >=0) counts[meta.unlim_idx] = 1;       // But if there's an unlimited dimension, set that count to 1
		strides.clear(); strides.resize(meta.dimnames.size(), 1); // set all strides to 1

	}

	void init_stream(NcStream& ncstream, std::string varname = ""){
		readMeta(ncstream.current_file, varname);
		current_file_id = ncstream.current_index.f_idx;
		streaming = true;
		if (t_idx < 0) throw std::runtime_error("GeoCube: Stream files must have a time dimension\n");
	}

	void print(bool b_values = false){
		std::cout << "Var: " << meta.name << " (" << meta.unit << ")\n";
		std::cout << "   dim names: " << meta.dimnames << '\n';
		std::cout << "   dim sizes (original): " << meta.dimsizes << '\n';
		std::cout << "      lat axis = " << lat_idx << "\n";
		std::cout << "      lon axis = " << lon_idx << "\n";
		std::cout << "      unlimited axis = ";
		if (meta.unlim_idx < 0) std::cout << "NA\n";
		else std::cout << meta.dimnames[meta.unlim_idx] << " (" << meta.unlim_idx << ")\n";
		std::cout << "   missing value = " << this->missing_value << "\n";
		std::cout << "   scale factor = " << scale_factor << "\n";
		std::cout << "   add offset = " << add_offset << "\n";
		std::cout << "   dimensions:\n";
		for (int i=0; i<meta.dimnames.size(); ++i){
			std::cout << "      " << meta.dimnames[i] << ": ";
			if (coords_trimmed[i].size() <= 6) std::cout << coords_trimmed[i] << '\n';
			else{
				auto& v = coords_trimmed[i];
				size_t n = coords_trimmed[i].size();
				std::cout << n << " | " << v[0] << " " << v[1] << " " << v[2] << " ... " 
				          << v[n-3] << " " << v[n-2] << " " << v[n-1] << "\n";
			}
		}

		Tensor<T>::print("   ", b_values);
	}

	// Note: Coord bounds are not updated when file is updated in streaming. This is
	// because a stream is supposed to have an identical grid across all files
	// Whether the grid is identical or not is not checked. It's up to the user to ensure this.
	void setCoordBounds(size_t axis, float lo, float hi){
		double descending = *coords[axis].rbegin() < *coords[axis].begin(); // check if the coordinate values are descending
		// get the start and end coordinate values that are just outside the lo-hi range
		size_t start=0, end=0;
		if (!descending){
			for (int i=0; i< coords[axis].size(); ++i) if (coords[axis][i] < lo) start = i; 
			for (int i=coords[axis].size()-1; i>=0; --i) if (coords[axis][i] > hi) end = i; 
		}
		else{
			for (int i=0; i< coords[axis].size(); ++i) if (coords[axis][i] > hi) start = i; 
			for (int i=coords[axis].size()-1; i>=0; --i) if (coords[axis][i] < lo) end = i; 
		}
		starts[axis] = start;
		counts[axis] = end - start + 1;

		coords_trimmed[axis].assign(coords[axis].begin()+start, coords[axis].begin()+end+1);

		// std::cout << "axis = " << dimnames[axis] << '\n';
		// std::cout << "Start: " << start << " " << coords[axis][start] << '\n';
		// std::cout << "End: " << end << " " << coords[axis][end] << '\n';
		// if (!descending) std::cout << "<<< " << coords[axis][start] << " | " << lo <<  " " << hi << " | " << coords[axis][end] << "\n";
		// else std::cout << ">>> " << coords[axis][start] << " | " << hi <<  " " << lo << " | " << coords[axis][end] << "\n";
	}

	void setIndices(size_t axis, size_t _start, size_t _count, ptrdiff_t _stride = 1){
		starts[axis]  = _start;
		counts[axis]  = _count;
		strides[axis] = _stride;
	}

	void readBlock(size_t unlim_start, size_t unlim_count){
		if (meta.unlim_idx >= 0){
			starts[meta.unlim_idx] = unlim_start;
			counts[meta.unlim_idx] = unlim_count;
		}
		std::cout << "Resizing tensor to: " << counts << '\n';
		this->resize(counts);
		ncvar.getVar(starts, counts, strides, this->vec.data());
	}

	void streamBlock(NcStream& ncstream, double julian_day){
		StreamIndex sid = ncstream.julian_to_indices(julian_day);
		// TODO: Check whether sid is same as current index, and if so, skip reading

		// if desired time is in not in current file, update file
		if (ncstream.current_index.f_idx != sid.f_idx){
			ncstream.update_file(sid.f_idx);
		}

		// if stream file has changed compared to the file last read in this->meta, update metadata
		// TODO: This can be put in above condition after update_file()
		if (current_file_id != sid.f_idx){
			init_stream(ncstream, meta.name);
		}

		starts[t_idx] = sid.t_idx;
		counts[t_idx] = 1;

		std::cout << "Resizing tensor to: " << counts << '\n';
		std::cout << "Reading frame @ " << ncstream.streamIdx_to_datestring(sid) << '\n';
		this->resize(counts);
		ncvar.getVar(starts, counts, strides, this->vec.data());
	}


};

} // namespace flare

#endif
