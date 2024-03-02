#ifndef MATH_TENSOR_H_
#define MATH_TENSOR_H_

#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>
#include <string>

/**
 Tensor. Multidimensional Array
 
 A tensor with dimensions
 ```
 { N2, N1, N0 } 
 ```
 is as follows.
 \image html tensorlib_tensor.png width=400cm 
 
 Elements of this tensor are accessed with a vector of indices as
 ```
 { i2, i1, i0 }
 ```
 i0 is the lowest dimension, i.e., (elements along i0 are stored consequtively in memory.
  This order of indices is chosen rather than {i0, i1, i2} to allow matrices 
  (two dimensional tensors) to be refered as {irow, icolumn}. 
*/
// Thus a 2x3x5 tensor will look like this:
// Tensor:
//   axis   2 1 0
//   dims = 2 3 5 
//   vals = 
//	  +-----+--+--+--+--+------> axis 0
//    |+    0  1  2  3  4 
//    | \   5  6  7  8  9 
//    |  \ 10 11 12 13 14 
//    v   \ .                        
//  axis 1 +     15 16 17 18 19 
//          \    20 21 22 23 24 
//           \   25 26 27 28 29 
//     axis 2 V
/**
 The indices of each element of a 2x3x5 Tensor are 
 ``` 
   loc: index 
   ---:------
     0: 0 0 0 
     1: 0 0 1 
     2: 0 0 2 
     3: 0 0 3 
     4: 0 0 4 
     5: 0 1 0 
     6: 0 1 1 
          :   
    25: 1 2 0 
    26: 1 2 1 
    27: 1 2 2 
    28: 1 2 3 
    29: 1 2 4  
 ```   
 
 */

// template <class T>
// std::ostream& operator << (std::ostream& os, const std::vector<T>& vec){
// 	for (auto v : vec) os << v << ' ';
// 	return os;
// }

template <class T>
class Tensor{
	private:
	std::vector<int> offsets;
	int nelem;
	
	public:
	std::vector<int> dim;
	std::vector<T> vec;

	// metadata
	T missing_value = std::numeric_limits<T>::quiet_NaN();

	public:

	Tensor<T>(){
		nelem = 0;
	}

	Tensor<T>(std::vector<int> _dim){
		resize(_dim);
	}

	const std::vector<int>& get_offsets() const {
		return offsets;
	}

	/// Create a tensor with specified dimensions.
	/// This function also allocates space for the tensor, and calculates the offsets used for indexing.
	void resize(std::vector<int> _dim){
		dim = _dim;
		nelem = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<int>());
		vec.resize(nelem);

		int ndim = dim.size();
		offsets.resize(ndim,0);
		int p = 1;
		for (int i=ndim-1; i>=0; --i){
			offsets[i] = p;
			p *= dim[i];
		}
		
	}

	void resize(std::vector<size_t> _dim){
		std::vector<int> dim1(_dim.begin(), _dim.end());
		resize(dim1);
	}

	/// Print the tensor.
	/// If vals is true, then values are also printed. Otherwise, only metadata is printed.
	void print(std::string prefix, bool vals){
	    std::cout << prefix << "Tensor:\n";
	    std::cout << prefix << "   dims = "; for (auto d : dim) std::cout << d << " "; std::cout << "\n";
	    std::cout << prefix << "   offs = "; for (auto d : offsets) std::cout << d << " "; std::cout << "\n";
	    std::cout << prefix << "   missing value = " << missing_value << "\n";
		if (vals){
			std::cout << prefix << "   vals = \n      "; std::cout.flush();
			for (int i=0; i<nelem; ++i){
				std::cout << prefix; 
				if (vec[i] == missing_value) std::cout << "NA";
				else if (std::isnan(missing_value) && std::isnan(vec[i])) std::cout << "NA";
				else std::cout << vec[i];
				std::cout << ' '; 
				bool flag = true;
				for (int axis=dim.size()-1; axis>0; --axis){
					flag = flag && (index(i)[axis] == dim[axis]-1);
					if (flag) std::cout << "\n      ";
				}
			}
		}
		std::cout << "\n";
	}

	void print(bool vals = true){
		print("", vals);
	}


//	TODO: 
//	This function can be private
	/// Convert coordinates (specified as a vector of indices) to 1D index where the value resides in the underlying vector.
	int location(std::vector <int> ix){
		int loc = 0;
		int ndim = dim.size();
		for (int i=ndim-1; i>=0; --i){
			loc += offsets[i]*ix[i];
		}
		return loc;
	}

	template<class... ARGS>
	int location(ARGS... ids){
		return location({ids...});
	}


	/// @brief Get the value at coordinates specified as a comma separated list. 
	///        The order of coordinates is \f$\{i_n, i_{n-1}, ..., i_0\}\f$
	template<class... ARGS>
	T& operator() (ARGS... ids){
		return vec[location({ids...})];
	}

	/// @brief Get the value at coordinates specified as an integer vector. 
	///        The order of coordinates is \f$\{i_n, i_{n-1}, ..., i_0\}\f$.
	T& operator() (std::vector<int> ix){
		return vec[location(ix)];
	}


	/// Convert 1D index to coordinates (Inverse of location())
	std::vector<int> index(int loc) const {
		int ndim = dim.size();
		std::vector<int> id(ndim);
		for (int i=ndim-1; i>=0; --i){
			int ix = loc % dim[i];
			loc = (loc-ix)/dim[i];
			id[i]=ix;
		}
		return id;
	}


	/// A utility function for testing purposes. Fills the tensor with incremental integers. 
	void fill_sequence(){
		for(size_t i=0; i<vec.size(); ++i) vec[i]=i;
	}

	
	/// @brief generate a list of 1D indices corresponding to all points on the hyperplane 
	/// perpendicular to 'axis' located at index 'k' on the axis. The axis is specified
	/// as the position (counted from the right) of the corresponding dimension in the 
	/// dimensions vector, i.e., between [0, n-1].
	//  [..., n2, n1, n0] = dimensions vector
	//  [...,  2,  1,  0] 
	//             ^ axis
	std::vector<int> plane(int axis, int k = 0) const {
		axis = dim.size()-1-axis;
		std::vector<int> locs;
		locs.reserve(nelem);
		for (int i=0; i<nelem; ++i){
			if (index(i)[axis] == 0) locs.push_back(i + k*offsets[axis]);
		}
		return locs;
	}

	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	void transform_dim(int loc, int axis, BinOp binary_op, std::vector<double> w){
		assert(w.size() == dim[dim.size()-1-axis]);
		
		axis = dim.size()-1-axis;
		int off = offsets[axis];
		
		for (int i=loc, count=0; count<dim[axis]; i+= off, ++count){
			if (vec[i] != missing_value) vec[i] = binary_op(vec[i], w[count]);	// this order is important, because the operator may not be commutative
		}
		
	}

	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	Tensor<T>& transform(int axis, BinOp binary_op, std::vector<double> w){
		std::vector<int> locs = plane(axis);
		for (int i=0; i<locs.size(); ++i){
			transform_dim(locs[i], axis, binary_op, w);
		}
		return *this;
	}
	
	
	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	/// @brief  Weighted accumulation of vector `v` at location `loc` along given axis: accumulate(w*v)
	/// @tparam BinOp 
	/// @param v0   Starting value over which to accumulate
	/// @param loc  The location at which `axis` is anchored - must be at index zero along the `axis` dimension 
	/// @param axis Dimension along which to accumulate - counted from the right
	/// @param binary_op Operator to use to successively accumulate values
	/// @param weights weights vector which multiplies the vector prior to accumulation
	/// @return Accumulated value
	template <class BinOp>
	T accumulate_dim(T v0, int loc, int axis, BinOp binary_op, std::vector<double> weights={}){
		assert(weights.size() == 0 || weights.size() == dim[dim.size()-1-axis]);
		
		axis = dim.size()-1-axis;
		int off = offsets[axis];
		
		T v = v0;
		for (int i=loc, count=0; count<dim[axis]; i+= off, ++count){
			double w = (weights.size()>0)? weights[count] : 1;
			if (vec[i] != missing_value) v = binary_op(v, w*vec[i]);
		}
		
		return v;
	}


	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	Tensor<T> accumulate(T v0, int axis, BinOp binary_op, std::vector<double> weights={}){
		std::vector<int> dim_new = dim;
		dim_new.erase(dim_new.begin()+dim_new.size()-1-axis);
		Tensor<T> tens(dim_new);
		tens.missing_value = this->missing_value;
		
		std::vector<int> locs = plane(axis);
		
		for (size_t i=0; i<locs.size(); ++i){
			tens.vec[i] = accumulate_dim(v0, locs[i], axis, binary_op, weights);
		}
		
		return tens;
	}


	Tensor<T> max_dim(int axis){
		T v0 = vec[1];
		return accumulate(v0, axis, [this](T a, T b){
			                            return (b == this->missing_value)? this->missing_value : std::max(a,b);
			                            });
	}


	int count_non_missing_dim(int loc, int axis){
		axis = dim.size()-1-axis;
		int off = offsets[axis];
		
		int n = 0;
		for (int i=loc, count=0; count<dim[axis]; i+= off, ++count){
			if (vec[i] != missing_value) ++n;
		}
		
		return n;
	}

	// count number of non-missing values along a given axis
	// useful for calculating averages along a dimension
	Tensor<int> count_non_missing(int axis){
		std::vector<int> dim_new = dim;
		dim_new.erase(dim_new.begin()+dim_new.size()-1-axis);
		Tensor<int> tens(dim_new);
		tens.missing_value = this->missing_value;
		
		std::vector<int> locs = plane(axis);
		
		for (size_t i=0; i<locs.size(); ++i){
			tens.vec[i] = count_non_missing_dim(locs[i], axis);
		}
		
		return tens;
	}

	Tensor<T> avg_dim(int axis){
		Tensor<T> tens = accumulate(0, axis, std::plus<T>());
		Tensor<int> counts = count_non_missing(axis);
		tens /= counts;
		return tens;
	}


	Tensor<T> repeat_inner(int n) const {
		std::vector<int> dim_new = dim;
		dim_new.push_back(n);
		
		Tensor<T> tout(dim_new);
		tout.missing_value = this->missing_value;
		int count = 0;
		for (int i=0; i<nelem; ++i){
			for (int j=0; j<n; ++j){
				tout.vec[count++] = vec[i];
			}
		}

		return tout;
	}

	Tensor<T> repeat_outer(int n) const {
		std::vector<int> dim_new = dim;
		dim_new.insert(dim_new.begin(), n);
		
		Tensor<T> tout(dim_new);
		tout.missing_value = this->missing_value;
		int count = 0;
		for (int j=0; j<n; ++j){
			for (int i=0; i<nelem; ++i){
				tout.vec[count++] = vec[i];
			}
		}

		return tout;
	}

	template <class S>
	Tensor<T> append(const Tensor<S> &rhs, int axis){
		int axis_right = axis;
		int axis_left = dim.size()-1-axis;
		std::vector<int> dim_new = this->dim;
		// check that dimensions conform along all axes except appending dimension
		for (int i=0; i<dim_new.size(); ++i){
			if (i == axis_left) continue;
			assert (dim_new[i] == rhs.dim[i]);
		}

		dim_new[axis_left] += rhs.dim[axis_left];
		Tensor<T> res(dim_new);
		res.missing_value = this->missing_value;

		int off_lhs = this->offsets[axis_left];
		int off_rhs = rhs.get_offsets()[axis_left];
		int off_res = res.get_offsets()[axis_left];

		// Note plane function requires axis counted from the right
		std::vector<int> locs_res1 = res.plane(axis_right); // get locations at index 0 on axis
		std::vector<int> locs_lhs = this->plane(axis_right);
		std::vector<int> locs_rhs = rhs.plane(axis_right);
		// for each location, firts copy values from lhs, then from rhs
		for (int il=0; il < locs_res1.size(); ++il){
			// at each location, copy lhs size worth of data 
			int i_res=locs_res1[il], i_lhs=locs_lhs[il], i_rhs=locs_rhs[il], count=0;
			for (; count<this->dim[axis_left]; i_res+= off_res, i_lhs += off_lhs, ++count){
				res.vec[i_res] = (this->vec[i_lhs] == this->missing_value)? res.missing_value : this->vec[i_lhs];	
			}
			// continue to copy rhs size worth of data 
			for (; count<dim_new[axis_left]; i_res+= off_res, i_rhs += off_rhs, ++count){
				res.vec[i_res] = (rhs.vec[i_rhs] == rhs.missing_value)? res.missing_value : rhs.vec[i_rhs];	
			}
		}

		return res;
	}


	Tensor<T> slice(int axis, size_t start, size_t end){
		int axis_right = axis;
		int axis_left = dim.size()-1-axis;
		std::vector<int> dim_new = this->dim;
		dim_new[axis_left] = end-start+1;

		Tensor<T> res(dim_new);
		res.missing_value = this->missing_value;

		int off_lhs = this->offsets[axis_left];
		int off_res = res.get_offsets()[axis_left];

		// Note plane function requires axis counted from the right
		std::vector<int> locs_res1 = res.plane(axis_right); // get locations at index 0 on axis
		std::vector<int> locs_lhs = this->plane(axis_right, start);
		// for each location, copy values from lhs
		for (int il=0; il < locs_res1.size(); ++il){
			// at each location, copy data from lhs
			int i_res=locs_res1[il], i_lhs=locs_lhs[il], count=0;
			for (; count < dim_new[axis_left]; i_res+= off_res, i_lhs += off_lhs, ++count){
				res.vec[i_res] = (this->vec[i_lhs] == this->missing_value)? res.missing_value : this->vec[i_lhs];	
			}
		}

		return res;
	}

	// Reverse along the given axis
	Tensor<T>& reverse(int axis){
		int axis_right = axis;
		int axis_left = dim.size()-1-axis;

		std::vector<int> locs = plane(axis_right);
		int off = offsets[axis_left];
		int N   = dim[axis_left];

		for (int il=0; il < locs.size(); ++il){
			for (int i=0; i<N/2; ++i){
				std::swap(vec[locs[il]+i*off], vec[locs[il]+(N-1-i)*off]);
				// std::cout << "swapping " <<  vec[locs[il]+i*off] << " and " << vec[locs[il]+(N-1-i)*off] << '\n';
			}
		}

		return *this;
	}	


	// Shift cyclically along the given axis by M indices
	// Right Shift by M (here = 4)
	//  1 2 3  4 | 5 6 7 8 9 10    <-- original vector
	// 10 9 8  7 | 6 5 4 3 2  1    <-- reverse entire range [0 --- N-1]
	//  7 8 9 10 | 1 2 3 4 5  6    <-- reverse ranges [0 --- M-1], [M --- N-1]
	Tensor<T>& rotate(int axis, int M){
		int axis_right = axis;
		int axis_left = dim.size()-1-axis;

		std::vector<int> locs = plane(axis_right);
		int off = offsets[axis_left];
		int N   = dim[axis_left];

		for (int il=0; il < locs.size(); ++il){
			// reverse [0 --- N-1]
			for (int i=0; i<N/2; ++i) std::swap(vec[locs[il]+i*off], vec[locs[il]+(N-1-i)*off]);
			// reverse [0 --- M-1]
			for (int i=0; i<M/2; ++i) std::swap(vec[locs[il]+i*off], vec[locs[il]+(M-1-i)*off]);
			// reverse [M --- N-1]
			for (int i=0; i<(N-M)/2; ++i) std::swap(vec[locs[il]+(M+i)*off], vec[locs[il]+(N-1-i)*off]);
		}

		return *this;
	}


	// set all values where msk is 0 or missing, to missing value
	template <class S>
	Tensor<T>& mask(const Tensor<S>& msk){
		assert(dim == msk.dim);

		auto mask_lambda = [this, &msk](T x, S m){
			return (m == 0 || m == msk.missing_value)? this->missing_value : x;
		};
		std::transform(vec.begin(), vec.end(), msk.vec.begin(), vec.begin(), mask_lambda);
		return *this;
	}

	// set all values where f returns false to missing value
	template <class Func>
	Tensor<T>& mask(Func f){
		auto mask_lambda = [this, &f](T x){
			return (f(x))? x : this->missing_value; // if function returns false, set to missing
		};
		std::transform(vec.begin(), vec.end(), vec.begin(), mask_lambda);
		return *this;
	}


	// Operators
	// see https://stackoverflow.com/questions/4421706/what-are-the-basic-rules-and-idioms-for-operator-overloading/4421719#4421719
	// All operators are non-commutative in the sense that the missing value and any other metadata is retained from the LHS
	public: 	
	template <class S>
	Tensor<T>& operator += (const Tensor<S>& rhs){
		assert(dim == rhs.dim);

		auto plus_missing = [this, &rhs](T x, S y){
			return (x == this->missing_value || y == rhs.missing_value)? this->missing_value : (x + y);
		};
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), plus_missing); //std::plus<T>());
		return *this;
	}
	
	template <class S>
	Tensor<T>& operator -= (const Tensor<S>& rhs){
		assert(dim == rhs.dim);

		auto minus_missing = [this, &rhs](T x, S y){
			return (x == this->missing_value || y == rhs.missing_value)? this->missing_value : (x - y);
		};
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), minus_missing); //std::minus<double>());
		return *this;
	}

	template <class S>
	Tensor<T>& operator *= (const Tensor<S>& rhs){
		assert(dim == rhs.dim);

		auto multiplies_missing = [this, &rhs](T x, S y){
			return (x == this->missing_value || y == rhs.missing_value)? this->missing_value : (x * y);
		};
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), multiplies_missing); //std::multiplies<T>());
		return *this;
	}

	/// @brief  Element wise division (this = this/rhs)
	/// @tparam S 
	/// @param rhs Tensor to divide by 
	template <class S>
	Tensor<T>& operator /= (const Tensor<S>& rhs){
		assert(dim == rhs.dim);

		auto divides_missing = [this, &rhs](T x, S y){
			return (x == this->missing_value || y == rhs.missing_value)? this->missing_value : (x / y);
		};
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), divides_missing);
		return *this;
	}

	template<class S>	
	Tensor<T>& operator += (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), 
		               [&s, this](const T& x){return (x == this->missing_value)? this->missing_value : (x+s);}
					   );
		return *this;
	}

	template <class S>
	Tensor<T>& operator -= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), 
		               [&s, this](const T& x){return (x == this->missing_value)? this->missing_value : (x-s);}
					   );
		return *this;
	}

	template <class S>
	Tensor<T>& operator *= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), 
		               [&s, this](const T& x){return (x == this->missing_value)? this->missing_value : (x*s);}
					   );
		return *this;
	}

	template<class S>
	Tensor<T>& operator /= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), 
		               [&s, this](const T& x){return (x == this->missing_value)? this->missing_value : (x/s);}
					   );
		return *this;
	}

};


template<class T, class S>
Tensor<T> operator + (Tensor<T> lhs, const Tensor<S>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs += rhs;
	return lhs;
}

template<class T, class S>
Tensor<T> operator - (Tensor<T> lhs, const Tensor<S>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs -= rhs;
	return lhs;
}

template<class T, class S>
Tensor<T> operator * (Tensor<T> lhs, const Tensor<S>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs *= rhs;
	return lhs;
}

template<class T, class S>
Tensor<T> operator / (Tensor<T> lhs, const Tensor<S>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs /= rhs;
	return lhs;
}

template<class T, class S>
Tensor<T> operator + (Tensor<T> lhs, S s){
	lhs += s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator - (Tensor<T> lhs, S s){
	lhs -= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator / (Tensor<T> lhs, S s){
	lhs /= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator * (Tensor<T> lhs, S s){
	lhs *= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator + (S s, Tensor<T> t){	// TODO: can passing be ref be used for these?
	return t+s;
}

template<class T, class S>
Tensor<T> operator - (S s, Tensor<T> t){
	return t-s;
}

template<class T, class S>
Tensor<T> operator * (S s, Tensor<T> t){
	return t*s;
}




#endif
