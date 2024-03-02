#ifndef FLARE_FLARE_CSVSTREAM_H
#define FLARE_FLARE_CSVSTREAM_H

#include <fstream>
#include <sstream>
#include <queue>

#include "csvrow.h"
#include "stream.h"

namespace flare{

class CsvStream : public Stream{
	protected:
	int t_id = -1;    ///< index at which time column exists in file
	bool store_data;  ///< Whether this stream stores data
	
	public:
	std::vector<std::string> colnames;  ///< Column names in csv file
	std::vector<CSVRow> data;           ///< data rows read from the csv files
	CSVRow current_row;                 ///< Last read row

	protected:
	inline CsvStream(bool _store_data) : store_data(_store_data){
	}

	public:
	// CsvStream constructor sets store_data to true
	inline CsvStream() : CsvStream(true){
	}

	inline void print_meta() override{
		Stream::print_meta();
		std::cout << "   colnames: " << colnames << '\n';
		std::cout << "   t_id: " << t_id << '\n';
	}

	inline void print_values(){
		std::cout << "time\t" << colnames << '\n';
		if (!store_data) return;
		for (int i=0; i<data.size(); ++i){
			std::cout << times[i] << '\t' << data[i] << '\n';
		}
	}

	inline void reset() override{
		Stream::reset();
		current_index.set(0,0,0);
		t_id = -1;
		colnames.clear();
		data.clear();
	}

	/// @brief Specializatin of Stream::open() for CSV files
	/// @param _filenames list of files to stream from
	inline void open(std::vector<std::string> _filenames, std::string _tunit_str){
		reset();

		filenames = _filenames;
		std::ifstream fin;
		
		// To construct the full time vector, we need to open each file once and obtain its time vector
		for (size_t i_file = 0; i_file<filenames.size(); ++i_file){
			std::string fname = filenames[i_file];
			fin.open(fname.c_str());
			if (!fin) throw std::runtime_error("Could not open file: "+fname);

			// read header
			CSVRow row;
			fin >> row;
			// store header only from 1st row of first file
			if (i_file == 0){
				for (int i=0; i<row.size(); ++i){
					colnames.push_back(row[i]);
				}
				// std::cout << colnames << std::endl;
			}

			// get index of time column (name comparisons are case insensitive)
			for (size_t icol=0; icol<colnames.size(); ++icol){
				std::string col = colnames[icol];
				for (auto name : tnames){
					if (utils::to_lower(col) == utils::to_lower(name)) t_id = icol;
				}
			}
			if (t_id < 0) throw std::runtime_error("Cannot find time column in CSV file: "+fname);

			idx_f0.push_back(times.size());
	
			// read time column for all rows
			int line_num = 0;
			while(fin >> row){
				times.push_back(std::stod(row[t_id]));
				if (store_data) data.emplace_back(row);
				file_indices.push_back(i_file);
				t_indices.push_back(line_num);
				++line_num;
			}

			fin.close();
		}

		// Parse time unit provided
		unit_str = _tunit_str;
		// std::cout << unit_str << std::endl;
		parse_time_unit(unit_str); // sets tunit, tscale, t_base

		if (!std::is_sorted(times.begin(), times.end()))  throw std::runtime_error("NcStream: Combined time vector is not in ascending order. Check the order of files supplied.\n");

		for (auto& t : times) t *= tscale; // convert time vector to "days since base date"

		if (times.size() > 0){
			tstep = (times.back() - times.front())/(times.size()-1); // get timestep in days
			DeltaT = times[times.size()-1] - times[0] + tstep;
		}

		if (store_data) current_row = data[0];
	}

	inline void advance_to_time(double j) override {
		// get index to read
		StreamIndex new_idx = julian_to_indices(j);
		if (debug) std::cout << "advance from " << current_index.f_idx << "." << current_index.t_idx << " --> " << new_idx.f_idx << "." << new_idx.t_idx << '\n';
		
		// Skip reading if new index is not different from current
		if (current_index == new_idx) return;

		current_row = data[new_idx.idx];
		current_index = new_idx;
	}
};


} // namespace flare


#endif
