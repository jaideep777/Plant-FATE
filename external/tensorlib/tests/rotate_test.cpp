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

	Tensor<double> u({2,7,10});
	u.fill_sequence();
	u.missing_value = -1;
	u.print();

	u.rotate(0, 4);
	u.print();

	u.rotate(1, 3);
	u.print();

	u.rotate(2, 1);
	u.print();

	vector<double> expected1 = { 
		116, 117, 118, 119, 110, 111, 112, 113, 114, 115, 
		126, 127, 128, 129, 120, 121, 122, 123, 124, 125, 
		136, 137, 138, 139, 130, 131, 132, 133, 134, 135, 
		76, 77, 78, 79, 70, 71, 72, 73, 74, 75, 
		86, 87, 88, 89, 80, 81, 82, 83, 84, 85, 
		96, 97, 98, 99, 90, 91, 92, 93, 94, 95, 
		106, 107, 108, 109, 100, 101, 102, 103, 104, 105, 

		46, 47, 48, 49, 40, 41, 42, 43, 44, 45, 
		56, 57, 58, 59, 50, 51, 52, 53, 54, 55, 
		66, 67, 68, 69, 60, 61, 62, 63, 64, 65, 
		6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 
		16, 17, 18, 19, 10, 11, 12, 13, 14, 15, 
		26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 
		36, 37, 38, 39, 30, 31, 32, 33, 34, 35
	};
	if (!equals(u.vec, expected1)) return 1;

	return 0;
}

