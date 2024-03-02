#include "../include/tensor.h"
#include <iostream>

using namespace std;

// compile:  g++ -Wall -Wextra -o 1 masked_addition_test.cpp 

bool equals(double x1, double x2, double tol=1e-6){
	return (abs(x1-x2) < tol);
}

template <class T>
bool equals(vector<T> x1, vector<T> x2, double tol=1e-6){
	if (x1.size()!= x2.size()) return false;
	for (size_t i=0; i<x1.size(); ++i) if (abs(x1[i]-x2[i]) > tol) return false;
	return true;
}

int main(){

	Tensor<double> u({2,3,5});
	u.fill_sequence();
	u.missing_value = -1;
	u(1,1,1) = -1;
	u(1,2,2) = -1;
	u(0,0,0) = -1;
	u(0,0,1) = -1;
	u.print();

	Tensor<double> v = u.accumulate(0, 1, std::plus<double>());
	v.print();

	vector<double> expected1 = 
		{15, 17, 21, 24, 27, 
         60, 42, 39, 69, 72
		 };
	if (!equals(v.vec, expected1)) return 1;

	u.transform(2, std::plus<double>(), {1,1});
	u.print();

	expected1 = {
		-1,-1,3,4,5,
		6,7,8,9,10,
		11,12,13,14,15,

		16,17,18,19,20,
		21,-1,23,24,25,
		26,27,-1,29,30
	};
	if (!equals(u.vec, expected1)) return 1;

	return 0;
}

