// CSV reader code adapted from the accepted answer to this post: 
// https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c

#ifndef FLARE_FLARE_CSVROW_H
#define FLARE_FLARE_CSVROW_H

#include <ios>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

namespace flare{

class CSVRow{
	private:
	std::string       m_line;
	std::vector<int>  m_data;

	public:
	inline std::string operator[](std::size_t index) const {
		std::string s(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
		s.erase(remove( s.begin(), s.end(), '\"' ),s.end()); // remove quotes from string
		return s;
	}

	inline std::size_t size() const {
		return m_data.size() - 1;
	}

	inline void readNextRow(std::istream& str){
		std::getline(str, m_line);
		m_line.erase(remove( m_line.begin(), m_line.end(), '\r' ), m_line.end()); // remove cariage return from string, if present

		m_data.clear();
		m_data.emplace_back(-1);
		std::string::size_type pos = 0;
		while((pos = m_line.find(',', pos)) != std::string::npos){
			m_data.emplace_back(pos);
			++pos;
		}
		// This checks for a trailing comma with no data after it.
		pos   = m_line.size();
		m_data.emplace_back(pos);
	}

	inline std::string get_line_raw(){
		return m_line;
	}
};

inline std::istream& operator>>(std::istream& str, flare::CSVRow& data){
    data.readNextRow(str);
    return str;
}   

inline std::ostream& operator<<(std::ostream& os, flare::CSVRow& data){
    for (size_t i=0; i<data.size(); ++i){
		os << data[i] << '\t';
	}
    return os;
}   

}

#endif
