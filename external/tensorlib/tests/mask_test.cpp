#include "../include/tensor.h"
#include <iostream>

using namespace std;

// compile:  g++ -Wall -Wextra -o 1 mask_test.cpp 

bool equals(double x1, double x2, double tol=1e-6){
	return (abs(x1-x2) < tol);
}

template <class T>
bool equals(vector<T> x1, vector<T> x2, double tol=1e-6){
	if (x1.size()!= x2.size()) return false;
	for (int i=0; i<x1.size(); ++i) if (abs(x1[i]-x2[i]) > tol) return false;
	return true;
}

int main(){

	{
		Tensor<double> u({2,3,5});
		u.fill_sequence();
		u.missing_value = -1;
		u.print();

		Tensor<double> msk({5});
		msk.vec = {1, 0, 0, 1, 0};
		msk.print();

		u.mask(msk.repeat_outer(3).repeat_outer(2));
		u.print();
	}

	{
		Tensor<double> u({2,3,5});
		u.fill_sequence();

		Tensor<double> msk({3});
		msk.vec = {1, 0, 0};
		msk.print();

		u.mask(msk.repeat_inner(5).repeat_outer(2));
		u.print();
	}


	return 0;
}
