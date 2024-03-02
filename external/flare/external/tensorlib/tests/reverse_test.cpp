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

	u.reverse(0);
	u.print();

	u.reverse(2);
	u.print();

	u.reverse(1);
	u.print();

	vector<double> expected1 = { 
		29, 28, 27, 26, 25, 
		24, 23, 22, 21, 20, 
		19, 18, 17, 16, 15, 

		14, 13, 12, 11, 10, 
		9, 8, 7, 6, 5, 
		4, 3, 2, 1, 0
	};
	if (!equals(u.vec, expected1)) return 1;

	return 0;
}

