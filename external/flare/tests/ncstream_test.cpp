#include "ncstream.h"
using namespace std;

int main(){

	flare::NcStream in_stream;
	in_stream.open({"tests/data/air.sig995.2011.nc",
	                "tests/data/air.sig995.2012.nc",
					"tests/data/air.sig995.2013.nc",
					"tests/data/air.sig995.2014.nc",
					"tests/data/air.sig995.2015.nc",
					});
	in_stream.print_meta();
	//in_stream.print_times();

	{
	/// TESTING OF T INDEX CALCULATION
	std::cout << "Calculate time index with periodic l-edged t vector\n\n";
	in_stream.periodic = true;
	in_stream.centered_t = false;
	double t0 = flare::datestring_to_julian("2010-12-30 00:00:00");
	for (double t = t0; t < t0+4; t+=0.125){
		std::cout << flare::julian_to_datestring(t) << ": " << in_stream.streamIdx_to_datestring(in_stream.julian_to_indices(t)) << "\n";
	}
	std::cout << "\n";

	std::cout << "Calculate time index with periodic l-edged t vector\n\n";
	t0 = flare::datestring_to_julian("2013-06-11 00:00:00");
	for (double t = t0; t < t0+2; t+=0.125){
		std::cout << flare::julian_to_datestring(t) << ": " << in_stream.streamIdx_to_datestring(in_stream.julian_to_indices(t)) << "\n";
	}
	std::cout << "\n";

	std::cout << "Calculate time index with periodic l-edged t vector\n\n";
	t0 = flare::datestring_to_julian("2015-12-30 00:00:00");
	for (double t = t0; t < t0+4; t+=0.125){
		std::cout << flare::julian_to_datestring(t) << ": " << in_stream.streamIdx_to_datestring(in_stream.julian_to_indices(t)) << "\n";
	}
	std::cout << "\n";
	}

	{
	/// TESTING OF T INDEX CALCULATION
	in_stream.periodic = true;
	in_stream.centered_t = true;
	std::cout << "Calculate time index with periodic centered t vector\n\n";
	double t0 = flare::datestring_to_julian("2010-12-30 00:00:00");
	for (double t = t0; t < t0+4; t+=0.125/4){
		std::cout << flare::julian_to_datestring(t) << ": " << in_stream.streamIdx_to_datestring(in_stream.julian_to_indices(t)) << "\n";
	}
	std::cout << "\n";
	}


	return 0;
}

