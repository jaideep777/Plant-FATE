#include "utils/moving_average.h"

using namespace std;

int main(){

	MovingAverager M;
	
	for (int i=0; i<10; ++i){
		M.push(i, 1);
		//M.print();
		cout << "t = " << i << ", Avg = " << M.get() << "\n";
	}

	M.clear();
	M.set_interval(1);
	
	for (int i=0; i<15; ++i){
		M.push(0.1*i, exp(-i/10.0));
		M.print();
		cout << "t = " << i << ", Avg = " << M.get() << "\n";
	}

	M.clear();
	M.set_interval(20);

	for (int i=0; i<100; ++i){
		M.push(i, (i<50)? 0.0 : 1.0);
		//M.print();
		cout << "t = " << i << ", Avg = " << M.get() << "\n";
	}


	return 0;
}

