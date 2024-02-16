#include "climate_stream.h"

namespace env{

void ClimateStream::init(){
	if (update_co2){
		co2_stream.set_tname("Year");
		co2_stream.periodic = false;
		co2_stream.centered_t = false;
		co2_stream.open({co2File}, "years CE");
		co2_stream.print_meta();
		co2_stream.print_times();
		co2_stream.print_values();
	}
	if (update_met){
		met_stream.set_tname("Decimal_year");
		met_stream.periodic = true;
		met_stream.centered_t = false;
		met_stream.open({metFile}, "years CE");
		met_stream.print_meta();
		met_stream.print_times();
		met_stream.print_values();
	}
}


void ClimateStream::updateClimate(double julian_time, Clim& C){
	if (update_co2){
		co2_stream.advance_to_time(julian_time);
		std::cout << co2_stream.current_row << std::endl;
		C.co2      = as<double>(co2_stream.current_row[1]); // co2 is in index 1
	}
	if (update_met){
		met_stream.advance_to_time(julian_time);
		std::cout << co2_stream.current_row << std::endl;
		C.tc       = as<double>(met_stream.current_row[3]);
		C.vpd      = as<double>(met_stream.current_row[4])*100; // convert hPa to Pa
		C.ppfd     = as<double>(met_stream.current_row[5]);
		C.ppfd_max = as<double>(met_stream.current_row[6]);
		C.swp      = as<double>(met_stream.current_row[7])*(-1); // convert -MPa to MPa
	}
}

}
