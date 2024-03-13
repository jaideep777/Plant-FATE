#ifndef FLARE_FLARE_TIMESTEPPER_H
#define FLARE_FLARE_TIMESTEPPER_H

#include <ios>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

#include "utils.h"
#include "time_math.h"

namespace flare{

/// @brief A universal stepper that can step in custom time units. 
/// Time units are specified as "units since base time". Here are a few examples:
/// 1. hours since 2020-01-01 0:0:0
/// 2. days since 1850-04-01 12:0:0
/// 3. months since 2000-01-01
/// 4. years since 1850-01-01
/// 3. years CE (equivalent to years since 0000-01-00)
class TimeStepper{
	public:
	double j_current;       ///< Current time [julian days]

	protected:
	std::string unit_str;   ///< full string representation of time unit (e.g., "days since yyyy-mm-dd hh:mm:ss")
	std::string tunit = ""; ///< time unit used by this stepper (e.g., "days", "months", etc)
	double tscale = 1;      ///< multiplier to convert time intervals from stepper's unit to 'days'
	std::tm t_base = {};    ///< epoch (base time) used by this stepper
	double j_base;          ///< epoch (base time) [julian days]

	public:
	inline void set_units(std::string time_units){
		parse_time_unit(time_units);
	}

	/// @brief     Advance current time by specified interval
	/// @param dt  time interval to advance, specified in stepper units
	/// @return   
	inline double step(double dt){
		double dt_days = dt*tscale; // convert specified interval into days
		j_current += dt_days;       // advance current julian date by specified days
		return j_current;
	}

	inline double to_julian(double t){
		return j_base + t*tscale;
	}

	inline double to_stepperUnits(double j){
		return (j - j_base)/tscale;
	}

	inline double get_tscale() const {
		return tscale;
	}

	protected:
	inline void parse_time_unit(std::string tunit_str){
		// treat "years CE" as "years since 0000-01-00 0:0:0". 
		// - This base date seems weird but works! 
		// - base date actually corresponds to -1-12-31 0:0:0, but negative years are not representable in the current format
		// - using "years since 0000-01-01" adds an extra day, 
		//       perhaps because over 2000 years (0001-2000) all leap days cancel out, 
		//       but 0000 is a leap year so adds an extra day
		if (tunit_str == "years CE") tunit_str = "years since 0000-01-00 0:0:0";

		// parse time units
		std::string since;
		std::stringstream ss(tunit_str);
		ss >> tunit >> since;

		if (since != "since") throw std::runtime_error("time unit is not in the correct format (<units> since <yyyy-mm-dd> <hh:mm:ss>)");

		if      (tunit == "days")     tscale = 1;
		else if (tunit == "hours")    tscale = 1.0/24.0;
		else if (tunit == "minutes")  tscale = 1.0/24.0/60.0;
		else if (tunit == "seconds")  tscale = 1.0/24.0/3600.0;
		else if (tunit == "months")   tscale = 1.0*365.2425/12;
		else if (tunit == "years")    tscale = 1.0*365.2425;

		if (tunit == "months" || tunit == "years") std::cout << "Warning: using " << tunit << " as time unit. 365.2425 days per year will be assumed. Conversion of time points to dates may have an error of +/- 1 day.\n";

		ss.str(tunit_str);
		ss >> std::get_time(&t_base, std::string(tunit + " since %Y-%m-%d %H:%M:%S").c_str());
		t_base.tm_zone = "GMT";

		j_base = date_to_julian(t_base);
	}

};

} // namespace flare

#endif
