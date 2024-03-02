#include "flare.h"

int main(){

	flare::NcStream in_stream;
	in_stream.periodic = true;
	in_stream.centered_t = true;

	in_stream.open({"tests/data/gpp.2000-2015.nc"});
	in_stream.current_file.printMeta();

	flare::GeoCube<float> v;
	v.init_stream(in_stream);
	v.print();

	v.setIndices(v.lat_idx, 230, 1);
	v.setIndices(v.lon_idx, 161, 2);
	v.readBlock(0,1);
	v.print(true);

	v.setCoordBounds(v.lat_idx, 18.5, 18.5);
	v.setCoordBounds(v.lon_idx, 77.25, 77.25);
	v.streamBlock(in_stream, flare::datestring_to_julian("2003-01-05 00:00:00"));
	v.print(true);



	return 0;
} 

