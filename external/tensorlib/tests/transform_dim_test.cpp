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

	u.transform(1, std::plus<double>(), {1,2,3});
	u.print();

	vector<double> expected_u = {
		1, 2, 3, 4, 5,      // +1
		7, 8, 9, 10, 11,    // +2
		13, 14, 15, 16, 17, // +3

		16, 17, 18, 19, 20, // +1
		22, 23, 24, 25, 26, // +2
		28, 29, 30, 31, 32  // +3
	};

	if (!equals(u.vec, expected_u)) return 1;

	u.transform(2, std::multiplies<double>(), {1,-1});
	u.print();

	vector<double> expected_u1 = {
		1, 2, 3, 4, 5,      // *1
		7, 8, 9, 10, 11,    // *1
		13, 14, 15, 16, 17, // *1

		-16, -17, -18, -19, -20, // *-1
		-22, -23, -24, -25, -26, // *-1
		-28, -29, -30, -31, -32  // *-1
	};

	if (!equals(u.vec, expected_u1)) return 1;

	return 0;
}

