#include <iostream>
#include <vector>
#include <string>
using namespace std;

using Matrix = std::vector<std::vector<string>>;
using Vector = std::vector<string>;

int main(){
	int K = 3;

	Matrix g_gx = {{"  g0", "  g1", "  g2"},
				   {"g0x0", "g1x0", "g2x0"},
				   {"g0x1", "g1x1", "g2x1"},
				   {"g0x2", "g1x2", "g2x2"},
					};

	Vector m_mx = {"   m", "m_x1", "m_x2", "m_x3"};

	cout << "g_gx:\n";
	for (int i=0; i<K+1; ++i){
		for (int j=0; j<K; ++j){
			cout << g_gx[i][j] << " ";
		}
		cout << '\n';
	}

	// This code goes into IEBT
	// ---------------------
	Matrix A(K+1);
	for (int i=0;i<K+1; ++i) A[i].resize(K+1, "----");

	for (int j=0; j<K; ++j){
		A[0][j+1] = m_mx[j+1]; // m_mx[1..K] is mortality gradient, m_mx[0] is the mortality rate
	}
	// all other rows (i=1..K), columns (j=0..K) contain the transposed g_gx matrix
	for (int i=0; i<K; ++i){ // row index goes from 0..K-1
		for (int j=0; j<K+1; ++j){ // column index goes from 0..K
			A[i+1][j] = g_gx[j][i]; 
			// ^ i+1 here because we need to skip the first row in A
		}
	}
	// ---------------------

	cout << "\nA:\n";
	for (int i=0; i<K+1; ++i){
		for (int j=0; j<K+1; ++j){
			cout << A[i][j] << " ";
		}
		cout << '\n';
	}

	Matrix expected = 
		{{"----", "m_x1", "m_x2", "m_x3"},
		 {"  g0", "g0x0", "g0x1", "g0x2"},
		 {"  g1", "g1x0", "g1x1", "g1x2"},
		 {"  g2", "g2x0", "g2x1", "g2x2"},
		};

	for (int i=0; i<K+1; ++i){
		for (int j=0; j<K+1; ++j){
			cout << A[i][j] << " " << expected[i][j] << "\n";
			if (A[i][j] != expected[i][j]) return 1;
		}
	}

	return 0;

}



