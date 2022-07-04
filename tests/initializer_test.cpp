#include "utils/initializer_v2.h"
//#include "assimilation.h"
//
using namespace std;

int main(){

	io::Initializer I;
	I.parse("tests/params/p_test_newsyntax.ini", false, true);

	I.print();

	cout << "~~~~~~~ retrieve values ~~~~~~~~~~~\n";
	cout << "sec1 / a = " << I.get<double>("sec 1", "a") << "\n";
	cout << "sec1 / b = " << I.get<double>("sec 1", "b") << "\n";
	cout << "sec2 / c = " << I.get<double>("sec2", "c") << "\n";

	vector<double> c = I.get_vector<double>("sec2", "c");
	cout << "sec2 / c = "; for (auto x : c) cout << x << ","; cout << "\n";

	vector<string> d = I.get_vector<string>("sec2", "d");
	cout << "sec2 / d = "; for (auto x : d) cout << x << ","; cout << "\n";

	vector<string> e = I.get_vector<string>("sec2", "e");
	cout << "sec2 / e = "; for (auto x : e) cout << x << ","; cout << "\n";

	cout << "glob / g1 = " << I.get<double>("g1") << "\n";
	cout << "glob / g2 = " << I.get<double>("global", "g2") << "\n";
	cout << "glob / g3 = " << I.get<string>("g3") << "\n";


	io::Initializer I2;
	stringstream sin;
	sin.str("# comment1\n\
	        [sec 1]\n\
	           a = 1  # value of a\n\
	           b= 2\n\
	        [sec2]\n\
	           #comment2\n\
	           ;comment3\n\
	           # comment 4\n\
	           c=3     8   2  1 4 1   6   ; haha\n\
	           d =4");
	I2.parse(sin, false, true);
	I2.print();

	return 0;
}

