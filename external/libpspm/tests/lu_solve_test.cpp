#include <iostream>
#include <cmath>
#include "../include/linear_system.h"
using namespace std;

// compile: g++ -o 1 lu_solve_test.cpp

int main() {

{
	Matrix A = {{1, 1, 0},
	            {2, 1, 3},
				{3, 1, 1}};

	Vector B = {1, 2, 5};			

	Vector X = luSolve(A,B);

	Vector expected_X =  {2.2, -1.2, -0.4};
	for (int i=0; i<X.size(); ++i){
		if (fabs(X[i] - expected_X[i])> 1e-6){
			std::cout << "FAIL\n";
			return 1;
		}
	}
}

{


	double m[3][3] = {{0, 1, 0},
	                  {2, 1, 3},
				      {3, 1, 1}};

	double * A[3] = {m[0], m[1], m[2]};

	cout << "double ** A:\n";
	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			cout << A[i][j] << ' ';
		}
		cout << '\n';
	}

	double B[3] = {1, 2, 5};			
	int P[4];
	

	LUPDecompose_wiki(A, 3, 1e-4, P);

	cout << "LU ** A:\n";
	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			cout << A[i][j] << ' ';
		}
		cout << '\n';
	}

	double X[3];
	LUPSolve_wiki(A, P, B, 3, X);

	cout << "X = "; 
	for (int i=0; i<3; ++i) std::cout << X[i] << ' ';
	cout << '\n';

	Vector expected_X =  {1.5714286,  1.0000000, -0.7142857};
	for (int i=0; i<3; ++i){
		if (fabs(X[i] - expected_X[i])> 1e-6){
			std::cout << "FAIL\n";
			return 1;
		}
	}

}


{
	Matrix A = {{0, 1, 0},
	            {2, 1, 3},
				{3, 1, 1}};

	Vector B = {1, 2, 5};			

	Vector X = lupSolve(A,B);

	Vector expected_X =  {1.5714286,  1.0000000, -0.7142857};
	for (int i=0; i<3; ++i){
		if (fabs(X[i] - expected_X[i])> 1e-6){
			std::cout << "FAIL\n";
			return 1;
		}
	}
}


}
