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
	u.print();

	Tensor<double> u0 = u.slice(0, 2,3);
	u0.print();

	vector<double> expected_u0 = {
		2, 3, 
		7, 8, 
		12, 13,
      
		17, 18, 
		22, 23, 
		-1., 28, 
	};

	if (!equals(u0.vec, expected_u0)) return 1;


	Tensor<double> u1 = u.slice(1, 0,1);
	u1.print();

	vector<double> expected_u1 = {
		0, 1, 2, 3, 4,       
		5, 6, 7, 8, 9,       
      
		15, 16, 17, 18, 19,  
		20, -1., 22, 23, 24, 
	};

	if (!equals(u1.vec, expected_u1)) return 1;

	Tensor<double> u2 = u.slice(2, 1, 1);
	u2.print();

	vector<double> expected_u2 = {
		15, 16, 17, 18, 19, 
		20, -1., 22, 23, 24, 
		25, 26, -1., 28, 29, 
	};

	if (!equals(u2.vec, expected_u2)) return 1;

	Tensor<double> u210 = u.slice(0, 1,2).slice(1, 0,1).slice(2, 0,0);
	u210.print();

	vector<double> expected_u210 = {
		1,2,
		6,7 
	};

	if (!equals(u210.vec, expected_u210)) return 1;

	return 0;
}

