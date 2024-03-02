#ifndef PSPM_FOF_H_
#define PSPM_FOF_H_

#include <vector>
#include <cmath>

#include "cohort.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Union-Find functions to calculate group Indices 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// given array of parents, find the root of q
inline int root(int q, std::vector<int>& par){
	while (q != par[q]){
		par[q] = par[par[q]];
		q = par[q];
	}
	return q;
}

// check if p and q have same root, i.e. belong to same group
inline bool find(int p, int q, std::vector<int>& par){
	return root(p, par) == root(q, par);
}

// put p and q in same group by merging the smaller tree into large one
inline void unite(int p, int q, std::vector<int>& par, std::vector<int>& sz){
	int i = root(p, par);
	int j = root(q, par);
	if (i==j) return;	// if both already have the same root do nothing
	if (sz[i] < sz[j]) {par[i]=j; sz[j] += sz[i];}
	else 			   {par[j]=i; sz[i] += sz[j];}
}


template<class Ind>
void group_cohorts(std::vector<Cohort<Ind>>& cohorts, double Rg){
	size_t N = cohorts.size();
	// std::cout << "N = " << N << '\n';
	std::vector <int> par(N);   // vector of parents
	std::vector <int> sz(N,1);  // vector of group sizes

	for (int i=0; i<N; ++i) par[i] = i;  // Set all parents to a unique number, i.e. all cohorts are in a different group
	
	long long int pairs = 0;  // number of pairs compared
	
	// FIXME: This O(N2) operation should be reduced by using spatial hashing or a kd-tree,
	//        but its not critical for O(100) cohorts
	for (int p=0; p<N; ++p){
		for (int q=0; q< p; ++q) {
			// std::cout << "try " << p << " and " << q << ": ";
			double d2 = 0;
			for (int k=0; k<cohorts[p].state_size; ++k){
				// FIXME: Maybe dx should this be normalized by maxSize()
				double dx = cohorts[p].x[k] - cohorts[q].x[k]; 
				d2 += dx*dx;
			}
			double d = sqrt(d2);
			// std::cout << "d = " << d << '\n';
			// if distance is < R_grp, assign same group
			if (d < Rg){
				// std::cout << "Uniting " << p << " and " << q << '\n';
				unite(p,q, par, sz);
			} 
			++pairs;
		}
	}

	for (int i=0; i<N; ++i){
		cohorts[i].group_id = root(i, par);
	}

	for (int i=0; i<N; ++i){
		cohorts[i].group_size = sz[cohorts[i].group_id];
	}

	// std::cout << "Pairs compared = " << pairs << std::endl;
	// std::cout << "Cohort group sizes: " << sz; 
	

}

#endif

