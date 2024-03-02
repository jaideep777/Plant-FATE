#ifndef PSPM_INDEX_UTILS_H_
#define PSPM_INDEX_UTILS_H_

#include <vector>
#include <numeric>

// From tensorlib

namespace utils{

namespace tensor{

/// Convert 1D index to coordinates (Inverse of location())
inline std::vector<int> index(int loc, const std::vector<int>& dim){
	int ndim = dim.size();
	std::vector<int> id(ndim);
	for (int i=ndim-1; i>=0; --i){
		int ix = loc % dim[i];
		loc = (loc-ix)/dim[i];
		id[i]=ix;
	}
	return id;
}


inline int location(const std::vector <int>& ix, const std::vector<int>& dim){
	int loc = 0;
	int ndim = dim.size();

	std::vector<int> offsets = std::vector<int>(ndim);
	int p = 1;
	for (int i=ndim-1; i>=0; --i){
		offsets[i] = p;
		p *= dim[i];
	}

	for (int i=ndim-1; i>=0; --i){
		loc += offsets[i]*ix[i];
	}
	return loc;
}

//       |
//  44   |
//  33   |------x  : id = [1,2]    
//  22   |      .    coords = [[10,20*,30,...]       coords[1]
//  11   |      .              [11,22,33*,44,...]]   coords[2]
//       |--------------------------------
//         10   20   30   40
inline std::vector<double> coord_value(const std::vector <int>& id, const std::vector<std::vector<double>>& coords){
	std::vector<double> vals(id.size());
	for (int i=0; i<id.size(); ++i){
		vals[i] = coords[i][id[i]];
	}
	return vals;
}

// check if index is a corner (n units inside the cube), i.e. if all index elements are n
inline bool is_corner_n(const std::vector <int>& ix, int n=0){
	return std::accumulate(ix.begin(), ix.end(), true, [n](const bool &acc, const int &j){return acc && (j==n);});
}

// check if index is on an (n units inside the cube), i.e. if any index element is n
inline bool is_edge_n(const std::vector <int>& ix, int n=0){
	return std::accumulate(ix.begin(), ix.end(), false, [n](const bool &acc, const int &j){return acc || (j==n);});
}

// Get the location of the index 1 less on axis k
inline int loc_minus1k(std::vector <int> ix, int k, const std::vector<int>& dim){
	ix[k] -= 1;
	assert(ix[k]>=0);
	return location(ix, dim);
}

// Get the location of the index 1 less on axis k
inline int loc_plus1k(std::vector <int> ix, int k, const std::vector<int>& dim){
	ix[k] += 1;
	assert(ix[k]<dim[k]);
	return location(ix, dim);
}


//       |
//  44   |      f14  f24
//  33   |------f13  f23     : id = [1,2]    
//  22   |      |            coords = [[10,20*,30,...]       coords[1]
//  11   |      |                     [11,22,33*,44,...]]   coords[2]
//  00   |--------------------------------
//         10   20   30   40
// 
//       id = [1,3]
//       g0 = (f23-f13)/h[0][1] = (f([2,3])-f([1,3])) / h[0][1]
//       g1 = (f14-f13)/h[1][2] = (f([1,4])-f([1,3])) / h[1][3]
// 
//       gk = (f([id+Ik]) - f([id])) / h[k][id[k]]
inline std::vector<double> grid_gradient(const std::vector <int>& id, const std::vector <int>& dim, const std::vector<double>& f, const std::vector<std::vector<double>>& x){
	std::vector<double> grad(id.size());
	for (int k=0; k<id.size(); ++k){
		std::vector<int>id_plus = id;
		id_plus[k] += 1;

		int j_plus   = location(id_plus, dim);
		int j_centre = location(id, dim);

		grad[k] = (f[j_plus] - f[j_centre]) / (x[k][id_plus[k]] - x[k][id[k]]);
	}
	return grad;
}

} // namespace tensor

namespace sequence{

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare){

	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	return p;
}

template <typename T>
void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p){
	std::vector<bool> done(vec.size());
	for (std::size_t i = 0; i < vec.size(); ++i){
		if (done[i]){
			continue;
		}
		done[i] = true;
		std::size_t prev_j = i;
		std::size_t j = p[i];
		while (i != j){
			std::swap(vec[prev_j], vec[j]);
			done[j] = true;
			prev_j = j;
			j = p[j];
		}
	}
}


inline std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

inline std::vector <double> logseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = exp(log(from) + i*(log(to)-log(from))/(len-1));
	return x;
}

inline std::vector <double> mids(const std::vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = (breaks[i]+breaks[i+1])/2;
	return a;
}

inline std::vector <double> left_edge(const std::vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = breaks[i];
	return a;
}

inline std::vector <double> right_edge(const std::vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = breaks[i+1];
	return a;
}

inline std::vector <double> diff(const std::vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = breaks[i+1]-breaks[i];
	return a;
}

} // namespace sequence

} // namespace utils


#endif
