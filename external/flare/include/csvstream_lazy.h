#ifndef FLARE_FLARE_CSVSTREAM_LAZY_H
#define FLARE_FLARE_CSVSTREAM_LAZY_H

#include <fstream>
#include <sstream>
#include <queue>

#include "csvrow.h"
#include "csvstream.h"

namespace flare{

class LazyCsvStream : public CsvStream{
	private:
	std::ifstream csvin;  ///< Input stream to read CSV file

	public:
	// LazyCsvStream constructor sets store_data to false
	inline LazyCsvStream() : CsvStream(false) {
	}

	/// @brief Specializatin of Stream::open() for CSV files
	/// @param _filenames list of files to stream from
	inline void open(std::vector<std::string> _filenames, std::string _tunit_str){
		CsvStream::open(_filenames, _tunit_str);
		// open first file
		update_file(0);

	}

	inline void update_file(size_t file_id){
		csvin.close();
		csvin.open(filenames[file_id]);
		if (!csvin) throw std::runtime_error("Could not open file: "+filenames[file_id]);

		csvin >> current_row; // skip header
		std::cout << " -- update_file(): header row: " << current_row.get_line_raw() << '\n';

		// set current index to first data entry in the file
		current_index.set(idx_f0[file_id], file_id, 0);
		std::cout << " -- updating file, index to --> " << current_index.f_idx << "." << current_index.t_idx << '\n';

		// Since data from current index must already be in current_row, read first data entry from file
		csvin >> current_row; 
	}


	inline void advance_to_time(double j) override {
		// get index to read
		StreamIndex new_idx = julian_to_indices(j);
		std::cout << "advance from " << current_index.f_idx << "." << current_index.t_idx << " --> " << new_idx.f_idx << "." << new_idx.t_idx << '\n';
		
		// Skip reading if new index is not different from current
		if (current_index == new_idx) return;

		// update file if 
		//    a) file index has changed
		//    b) new idx < current idx (since ifstream cannot go backwards)
		if (current_index.f_idx != new_idx.f_idx ||
		    current_index.idx > new_idx.idx){
			update_file(new_idx.f_idx);
		}

		std::cout << " -- advancing: " << current_index.f_idx << "." << current_index.t_idx << " --> " << new_idx.f_idx << "." << new_idx.t_idx << '\n';
		// consume rows until we hit new_index
		for (int i=current_index.t_idx; i<new_idx.t_idx; ++i){
			csvin >> current_row;
			std::cout << "   -- " << current_row.get_line_raw() << "\n";
		}
		current_index = new_idx;
	}
};


} // namespace flare


#endif
