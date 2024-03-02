#ifndef  PSPM_UTIL_LINEAR_SYSTEM_H_
#define  PSPM_UTIL_LINEAR_SYSTEM_H_

// Simple solver for linear system AX=B using LU decomposition
// Ref: https://en.wikipedia.org/wiki/LU_decomposition

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

inline void printMatrix(const Matrix& m, int n){
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
			std::cout<<m[i][j]<<"  ";
		std::cout<< '\n';
	}
	std::cout<< '\n';
}

// Code for this function is adapted from 
// http://blog.gtiwari333.com/2009/12/c-c-code-lu-decomposition-for-solving.html
inline Vector luSolve(const Matrix& a, const Vector& b){

	int n = a.size();
	assert(n == b.size());
	for (auto& v : a) assert(n == v.size());
	printMatrix(a,n);

	Matrix lu(n, Vector(n,0)); 
	printMatrix(lu,n);

	Vector z(n,0), x(n,0);

	//********** LU decomposition *****//
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[i][k] * lu[k][j];
			lu[i][j] = a[i][j] - sum;
		}
		for (int j = i + 1; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[j][k] * lu[k][i];
			lu[j][i] = (1 / lu[i][i]) * (a[j][i] - sum);
		}
	}
	//******** Displaying LU matrix**********//
	std::cout << "\nLU matrix is\n";
	printMatrix(lu,n);

	//***** Solve LZ=b*********//
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += lu[i][k] * z[k];
		z[i] = b[i] - sum;
	}    
	//***** Solve UX=Z*********//
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int k = i + 1; k < n; k++)
			sum += lu[i][k] * x[k];
		x[i] = (1 / lu[i][i]) * (z[i] - sum);
	}
	//*********** DISPLAYING SOLUTION**************//
	std::cout<<"\nSet of solution is\n";
	for(int i=0;i<n;i++)
		std::cout<<x[i]<<'\n';

	return x;
}


// Code for this function is from 
// https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
 */
inline int LUPDecompose_wiki(double **A, int N, double Tol, int *P) {

	int i, j, k, imax; 
	double maxA, *ptr, absA;

	for (i = 0; i <= N; i++)
		P[i] = i; //Unit permutation matrix, P[N] initialized with N

	for (i = 0; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) { 
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) return 0; //failure, matrix is degenerate

		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			ptr = A[i];
			A[i] = A[imax];
			A[imax] = ptr;

			//counting pivots starting from N (for determinant)
			P[N]++;
		}

		for (j = i + 1; j < N; j++) {
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}

	return 1;  //decomposition done 
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
inline void LUPSolve_wiki(double **A, int *P, double *b, int N, double *x) {

	for (int i = 0; i < N; i++) {
		x[i] = b[P[i]];

		for (int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int k = i + 1; k < N; k++)
			x[i] -= A[i][k] * x[k];

		x[i] /= A[i][i];
	}
}


// Code adapted from https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
inline Vector lupSolve(Matrix& A, const Vector& b, double Tol = 1e-6){
	int n = A.size();
	assert(n == b.size());
	for (auto& v : A) assert(n == v.size());
	// printMatrix(A,n);

	std::vector<int> P(n+1);
	std::iota(P.begin(), P.end(), 0);

	// Find LU
	for (int i = 0; i < n; i++) {
		double maxA = 0.0, absA;
		int imax = i;

		for (int k = i; k < n; k++)
			if ((absA = fabs(A[k][i])) > maxA) { 
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) throw std::runtime_error("failure, matrix is degenerate");

		if (imax != i) {
			//pivoting P
 			std::swap(P[i], P[imax]); 

			//pivoting rows of A
 			std::swap(A[i], A[imax]); // Swapping vectors is still O(1)

			//counting pivots starting from N (for determinant)
			P[n]++;
		}

		for (int j = i + 1; j < n; j++) {
			A[j][i] /= A[i][i];

			for (int k = i + 1; k < n; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}	

	// printMatrix(A, n);

	// Solve LUX = B
	Vector x(n);
	for (int i = 0; i < n; i++) {
		x[i] = b[P[i]];

		for (int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (int i = n - 1; i >= 0; i--) {
		for (int k = i + 1; k < n; k++)
			x[i] -= A[i][k] * x[k];

		x[i] /= A[i][i];
	}

	return x;
}

#endif