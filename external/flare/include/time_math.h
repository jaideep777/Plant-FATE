#ifndef FLARE_FLARE_TIME_MATH_H
#define FLARE_FLARE_TIME_MATH_H

#include <iomanip>
#include <iostream>
#include <chrono>
#include <sstream>

// Note: 
// std::tm format:
// year = years OVER 1900. So add 1900 to get CE year
// mon = 0-11. So add 1 for gday calculations
// mday = 1-31

// Reference: https://xkcd.com/2867/

namespace flare{

struct time_point {
	int year = 0;   // CE year = years over 0000
	int mon = 0;    // 1-12
	int day = 0;    // 1-31
	int hour = 0;   // 0-23
	int min = 0;    // 0-59
	int sec = 0;    // 0-59
};

inline bool isLeapYear(int year){
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

// time is assumed to be in GMT. Time zone conversion is up to the user
inline time_point string_to_date(const std::string& s, std::string format = "%d-%d-%d %d:%d:%d") {
	time_point t = {};

	int n = std::sscanf(s.c_str(), format.c_str(),
		&t.year, &t.mon, &t.day, &t.hour, &t.min, &t.sec);

	if (n < 3) throw std::invalid_argument("Invalid date string format");

	if (n < 6) t.hour = t.min = t.sec = 0;

	return t;
}

inline std::string date_to_string(const time_point& t) {
	char buf[64];
	std::sprintf(buf, "%04d-%02d-%02d %02d:%02d:%02d GMT",
			t.year, t.mon, t.day, t.hour, t.min, t.sec);
	return buf;
}

/// @brief calculate days since 01-Mar-0000 AD 
/// @param y year [anything]
/// @param m month [1-12]
/// @param d day [1-31]
/// @return days since 01-Mar-0000 AD 
/// see: http://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
/// archived here: https://web.archive.org/web/20170507133619/https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html 
inline int _ymd2gday(int y, int m, int d){
	m = (m+9)%12;
	y = y - m/10;
	int days_since_01_03_0000 = 365*y + y/4 - y/100 + y/400 + (m*306 + 5)/10 + (d-1);
	return days_since_01_03_0000; 
}


/// @brief Converts date struct to Julian day
/// @param time_struct time to convert
/// @return Julian day (including day fraction)
/// Julian days start on -4173-11-24 12:00:00 (24 Nov 4174 BC)
/// https://en.wikipedia.org/wiki/Julian_day
/// time is assumed to be in GMT. Time zone conversion is up to the user
inline double date_to_julian(time_point time_struct){
	double d_int = _ymd2gday(time_struct.year, time_struct.mon, time_struct.day);
	double d_frac = (double(time_struct.hour) + double(time_struct.min)/60.0 + double(time_struct.sec)/3600.0)/24.0;
	// std::cout << "d = " << d_int << " . " << d_frac << '\n';
	return d_int + d_frac + 1721119.5; // 1721119.5 is julian day number at 00:00:00 on 01-03-0000, i.e. the difference between gday epoch (01-03-0000 00:00:00) and julian epoch (-4173-11-24 12:00:00)
}


/// @brief Convert Julian day to time struct
/// @param julian_date Julian day (including day fraction)
/// @return time
/// First subtracts offset to convert from Julian day to 
///   days since 01-Mar-0000 AD (arbitrary reference date)
///   see: http://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
///   archived here: https://web.archive.org/web/20170507133619/https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html 
/// Then converts that to time struct
/// time is assumed to be in GMT. Time zone conversion is up to the user
inline time_point julian_to_date(double julian_date){
	time_point result;
	double gday = julian_date - 1721119.5;
	int g = int(gday);
	double dayf = gday - g;	// use of double here gives 5:29:59.9 for 5:30:0!!

	// get day in yyyy, mm, dd
	int ystr, mstr, dstr;
	int y, ddd, mm, mi;
	y = (10000*static_cast<long long int>(g) + 14780)/3652425;
	ddd = g - (365*y + y/4 - y/100 + y/400);
	if (ddd < 0){
		--y;
		ddd = g - (365*y + y/4 - y/100 + y/400);
	}
	mi = (52 + 100*ddd)/3060;
	result.year = ystr = y + (mi + 2)/12;
	result.mon = mstr = (mi + 2)%12 + 1;
	result.day = dstr = ddd - (mi*306 + 5)/10 + 1;
	
	// get time in hh, mm, ss
	dayf = dayf*24;
	result.hour = int(dayf);
	double r = dayf - int(dayf);
	result.min = int(r*60);
	r = r*60 - int(r*60);
	result.sec = r*60;

	return result;
}

inline std::string julian_to_datestring(double j){
	return date_to_string(julian_to_date(j));
}

inline double datestring_to_julian(std::string datestring){
	return date_to_julian(string_to_date(datestring));
}

inline double decimal_year(time_point tp){
	double days_in_yr = isLeapYear(tp.year)? 366:365;

	static const int daysToMonth[2][12] = {
		{ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
		{ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 },
	};
	
	double yday = daysToMonth[isLeapYear(tp.year)? 1 : 0][tp.mon-1] + tp.day-1;

	return tp.year + yday/days_in_yr;
}

inline double julian_to_yearsCE(double j){
	return (j - datestring_to_julian("0000-01-00 0:0:0"))/365.2425; 
}

inline double yearsCE_to_julian(double years_ce){
	return datestring_to_julian("0000-01-00 0:0:0") + years_ce*365.2425; 
}


} // namespace flare


#endif
