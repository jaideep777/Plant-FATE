#include <iostream>
#include <iomanip>
#include <pn_zero.h>
using namespace std;

int main(){
	auto f = [](double x){
		return pow((x-3),3) + (x-3)*(x-2);
	};
	
	auto res = pn::zero(2,8,f,1e-6);
	cout << "Root is: " <<  std::setprecision(20) << res.root << ", iter = " << res.nfnct <<  "\n";
	
	if (fabs(res.root - 3) > 1e-6) return 1;
	else return 0;
}
