#include "../include/tensor.h"
#include <iostream>

using namespace std;

// compile:  g++ -Wall -Wextra -o 1 masked_addition_test.cpp 

inline bool equals(double x1, double x2, double tol=1e-6){
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

	Tensor<double> v = u.accumulate(0, 0, std::plus<double>());
	vector<double> expected1;
	expected1 = {10,35,60,  
	            85,110,135};
	v.print();
	if (!equals(v.vec, expected1)) return 1;
	
	Tensor<double> v1 = u.accumulate(0, 1, std::plus<double>());
	expected1 = {15, 18, 21, 24, 27, 
                60, 63, 66, 69, 72};
	v1.print();
	if (!equals(v1.vec, expected1)) return 1;

	Tensor<double> v2 = u.accumulate(0, 2, std::plus<double>());
	expected1 = {15, 17, 19, 21, 23, 
                25, 27, 29, 31, 33, 
                35, 37, 39, 41, 43};
	v2.print();
	if (!equals(v2.vec, expected1)) return 1;


	return 0;
}

