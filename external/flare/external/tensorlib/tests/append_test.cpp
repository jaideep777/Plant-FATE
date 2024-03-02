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

	Tensor<int> v({2,2,5});
	v.fill_sequence();
	v += 100;
	v.missing_value = -2;
	v(0,0,0) = -2;
	v(0,0,1) = -2;
	v.print();

	Tensor<double> uv = u.append(v, 1);
	uv.print();

	vector<double> expected_uv = {
		0, 1, 2, 3, 4, 
		5, 6, 7, 8, 9, 
		10, 11, 12, 13, 14, 
		-1., -1., 102, 103, 104, 
		105, 106, 107, 108, 109, 
      
		15, 16, 17, 18, 19, 
		20, -1., 22, 23, 24, 
		25, 26, -1., 28, 29, 
		110, 111, 112, 113, 114, 
		115, 116, 117, 118, 119 
	};

	if (!equals(uv.vec, expected_uv)) return 1;

	Tensor<float> w({2,3,2});
	w.fill_sequence();
	w += 100;
	w.missing_value = -3;
	w(0,0,0) = -3;
	w(1,0,0) = -3;
	w.print();

	Tensor<double> uw = u.append(w, 0);
	uw.print();

	vector<double> expected_uw = {
		0, 1, 2, 3, 4,       -1., 101,
		5, 6, 7, 8, 9,       102, 103,
		10, 11, 12, 13, 14,  104, 105,
      
		15, 16, 17, 18, 19,   -1., 107,
		20, -1., 22, 23, 24,  108, 109,
		25, 26, -1., 28, 29,  110, 111
	};

	if (!equals(uw.vec, expected_uw)) return 1;

	Tensor<float> q({3,5});
	q.fill_sequence();
	q += 200;
	q.missing_value = -4;
	q(1,4) = -4;
	q(2,4) = -4;
	q.print();

	Tensor<double> uq = u.append(q.repeat_outer(1), 2);
	uq.print();

	vector<double> expected_uq = {
		0, 1, 2, 3, 4, 
		5, 6, 7, 8, 9, 
		10, 11, 12, 13, 14, 
      
		15, 16, 17, 18, 19, 
		20, -1., 22, 23, 24, 
		25, 26, -1., 28, 29, 

		200, 201, 202, 203, 204, 
		205, 206, 207, 208, -1.0, 
		210, 211, 212, 213, -1.0
	};

	if (!equals(uq.vec, expected_uq)) return 1;

	// u.append(v, 0); // should throw

	return 0;
}

