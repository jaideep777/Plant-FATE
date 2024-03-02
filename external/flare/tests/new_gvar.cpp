#include "flare.h"

int main(){

	flare::NcFilePP in_file;
	in_file.open("tests/data/gpp.2000-2015.nc", netCDF::NcFile::read);
	in_file.readMeta();
	in_file.printMeta();

	flare::GeoCube<float> v;
	v.readMeta(in_file);
	v.print();

	v.setIndices(v.lat_idx, 230, 1);
	v.setIndices(v.lon_idx, 161, 2);
	v.readBlock(0,1);
	v.print(true);

	// just a test of coord bounding
	std::reverse(v.coords[v.lat_idx].begin(), v.coords[v.lat_idx].end());
	v.setCoordBounds(v.lon_idx, 75, 101);
	v.setCoordBounds(v.lat_idx, 60, 80);

	v.setCoordBounds(v.lon_idx, 75.25, 101.75);
	v.setCoordBounds(v.lat_idx, 60.25, 80.75);

	v.setCoordBounds(v.lon_idx, 111.25, 131.75);
	v.setCoordBounds(v.lat_idx, -60.75, 12.25);
	std::reverse(v.coords[v.lat_idx].begin(), v.coords[v.lat_idx].end());
	// ~~~

	v.setCoordBounds(v.lat_idx, 18.5, 18.5);
	v.setCoordBounds(v.lon_idx, 77.25, 77.25);
	v.readBlock(0, 2);
	v.print(true);



	return 0;
} 

