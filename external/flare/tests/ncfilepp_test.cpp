#include "ncfilepp.h"
using namespace std;

int main(){

{
	flare::NcFilePP in_file;
	in_file.open("tests/data/gpp.2000-2015.nc", netCDF::NcFile::read);
	in_file.readMeta();
	in_file.printMeta();
}

{
	flare::NcFilePP in_file;
	in_file.open("tests/data/AR-SLu_2009-2011_FLUXNET2015_Met.nc", netCDF::NcFile::read);
	in_file.readMeta();
	in_file.printMeta();
}

{
	flare::NcFilePP in_file;
	in_file.open("tests/data/surta_india_0.2.nc", netCDF::NcFile::read);
	in_file.readMeta();
	in_file.printMeta();
}

	return 0;
}

