#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <netcdf>
#include <chrono>
#include <cmath>
#include "../include/time_math.h"
using namespace std;

// checks if 2 dates are within 1 sec of each other
bool second_equal(const std::tm& t1, const std::tm& t2){
	return 
		(t1.tm_hour == t2.tm_hour) &&
		(t1.tm_mday == t2.tm_mday) &&
		(t1.tm_min  == t2.tm_min)  &&
		(t1.tm_mon	== t2.tm_mon)  &&
		(t1.tm_year == t2.tm_year) &&
		(t1.tm_yday == t2.tm_yday) &&
		(fabs(t1.tm_sec  - t2.tm_sec) < 2);
}

int main(){
	vector<string> datestrings = {
		"2013-01-01 00:30:00", // https://en.wikipedia.org/wiki/Julian_day
		"2001-12-18 12:10:45",
		"1990-06-17 21:14:45",
		"1988-10-18 08:58:30", 
		"2988-07-01 12:00:00",
		"0000-03-01 12:00:00",
		"0000-03-01 00:00:00"
		// "0000-01-00 00:00:00"
		// "-2000-03-01 00:00:00",
		// "-4172-11-24 12:00:00"
	};

	// https://www.aavso.org/cgi-bin/jd2cal.pl
	vector<double> expected_julian = {
		2456293.52083, // https://en.wikipedia.org/wiki/Julian_day
		2452262.00747,
		2448060.38524,
		2447452.87396,
		2812587.00000,
		1721120.00000,
		1721119.50000
		// 1721059.50000
	};

	vector<std::tm> dates;
	for (auto s : datestrings){
		dates.push_back(flare::string_to_date(s));
		cout << flare::date_to_string(*dates.rbegin());
		// cout << " ("; print_date(*dates.rbegin()); cout << ")";
		cout << "\n";
	}

	vector<double> julians;
	for (auto d : dates){
		julians.push_back(flare::date_to_julian(d));
		cout << flare::date_to_string(d);
		cout << setprecision(12) << " = " << flare::date_to_julian(d) << '\n';
	}

	// Check date to julian conversion
	for (int i=0; i<julians.size(); ++i){
		cout << julians[i] << " === " << expected_julian[i] 
		     << " (diff = " << julians[i] - expected_julian[i] << ")\n";
		if (fabs(julians[i] - expected_julian[i]) > 1e-5){
			cout << "FAILED\n";
			return 1;
		}
	}

	vector<std::tm> dates_reconv;
	for (auto j : julians){
		cout << j << " = ";
		dates_reconv.push_back(flare::julian_to_date(j));
		cout << flare::date_to_string(*dates_reconv.rbegin());
		cout << '\n';
	}

	// Check julian to date conversion
	for (int i=0; i<julians.size(); ++i){
		cout << flare::date_to_string(dates[i]) << " === " << flare::date_to_string(dates_reconv[i]) << "\n";
		if (!second_equal(dates[i], dates_reconv[i])){
			cout << "FAILED\n";
			return 1;
		}
	}

	cout << "Julian @ 2000.00 = " << setprecision(12) << flare::yearsCE_to_julian(2000) << '\n';
	cout << "years CE @ 2000.00 = " << flare::julian_to_yearsCE(flare::yearsCE_to_julian(2000)) << '\n';
	if (fabs(flare::yearsCE_to_julian(2000) - 2451544.5) > 1e-6) return 1;
	if (fabs(flare::julian_to_yearsCE(flare::yearsCE_to_julian(2000))- 2000)>1e-6) return 1;

	cout << "-----------------\n";
	cout << "All tests PASSED!\n";

	return 0;
}

