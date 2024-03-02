#include "csvstream.h"
#include "csvstream_lazy.h"
using namespace std;

// to verify, run tests/Rscripts/csvstream_periodic_extension_test_analysis.R

template<class STREAM>
int test(string outfile){

	flare::CsvStream in_stream;
	in_stream.periodic = true;
	in_stream.centered_t = false;

	in_stream.open({"tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv"}, 
					"years CE");
	in_stream.print_meta();
	in_stream.print_times();

	in_stream.advance_to_time(flare::datestring_to_julian("2002-01-04"));
	in_stream.print_meta();
	cout << "Data:\n";
	// cout << '|' << in_stream.current_row.get_line_raw() << '|' << endl;
	for (int i=0; i<in_stream.current_row.size(); ++i) cout << in_stream.current_row[i] << "\t";
	cout << '\n';
	
	ofstream fout(outfile.c_str());
	double t0 = flare::datestring_to_julian("1921-01-04");
	double tf = flare::datestring_to_julian("2081-12-31");
	for (double t = t0; t <= tf; t += 365.2425/12.0/1.0){
		in_stream.advance_to_time(t);
		fout << t << "\t" << flare::decimal_year(flare::julian_to_date(t)) << "\t";
		for (int i=0; i<in_stream.current_row.size(); ++i) fout << in_stream.current_row[i] << "\t";
		fout << '\n';
	}
	fout.close();

	
	return 0;
}

int main(){

	cout << "-------------- CSV Stream ---------------\n\n";
	int res1 = test<flare::CsvStream>("csvstream_met.txt");

	cout << "----------- Lazy CSV Stream ---------------\n\n";
	int res2 = test<flare::LazyCsvStream>("csvstream_lazy_met.txt");

	return res1+res2;	
}

