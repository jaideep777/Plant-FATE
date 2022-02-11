#include "utils/initializer_v2.h"
//#include "assimilation.h"
//
using namespace std;

int main(){
	
	ifstream fin("tests/params/p_newsyntax.ini");
	io::Initializer I;
	I.parse(fin);

	return 0;
}

