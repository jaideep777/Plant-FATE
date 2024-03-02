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

	Tensor<float> v({2,3,5});
	v.fill_sequence();
	v.missing_value = -2;
	v(0,0,0) = -2;
	v(0,0,1) = -2;
	v.print();

	Tensor<double> sum = u + v;
	sum.print();

	vector<double> expected1 = 
		{-1, -1, 4, 6, 8, 
		 10, 12, 14, 16, 18, 
		 20, 22, 24, 26, 28,

		 30, 32, 34, 36, 38,
		 40, -1, 44, 46, 48,
		 50, 52, -1, 56, 58 
		 };
	if (!equals(sum.vec, expected1)) return 1;

	Tensor<double> diff = u - v;
	diff.print();

	expected1 = 
		{-1, -1, 0,0,0, 
		 0,0,0,0,0, 
		 0,0,0,0,0,

		 0,0,0,0,0, 
		 0,-1,0,0,0,
		 0,0,-1,0,0,
		 };
	if (!equals(diff.vec, expected1)) return 1;

	Tensor<double> div = (2*u) / v;
	div.print();

	expected1 = 
		{-1, -1, 2,2,2, 
		 2,2,2,2,2,
		 2,2,2,2,2,

		 2,2,2,2,2, 
		 2,-1,2,2,2,
		 2,2,-1,2,2
		 };
	if (!equals(div.vec, expected1)) return 1;


	return 0;
}

