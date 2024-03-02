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
	u.print();

	u.count_non_missing(2).print();

	u.max_dim(1).print();

	Tensor<double> v = u.avg_dim(2);
	v.print();

	vector<double> expected1 = { 
	   7.5,  8.5,  9.5, 10.5, 11.5, 
      12.5, 13.5, 14.5, 15.5, 16.5, 
      17.5, 18.5, 19.5, 20.5, 21.5
	};
	if (!equals(v.vec, expected1)) return 1;

	u(1,1,1) = -1;
	u(1,2,2) = -1;
	u(0,0,0) = -1;
	u(0,0,1) = -1;
	u.print();

	Tensor<double> umax = u.max_dim(1);
	umax.print();
	expected1 = {
		10, 11, 12, 13, 14, 
        25, 26, 22, 28, 29
	};
	if (!equals(umax.vec, expected1)) return 1;


	Tensor<double> v1 = u.avg_dim(2);
	v1.print();

	expected1 = { 
	    15,  16,  9.5, 10.5, 11.5, 
      12.5, 6, 14.5, 15.5, 16.5, 
      17.5, 18.5, 12, 20.5, 21.5
	};
	if (!equals(v1.vec, expected1)) return 1;


	// if (!equals(u.vec, expected1)) return 1;

	return 0;
}

