#ifndef ADAPTIVE_INTERPOLATOR_H_
#define ADAPTIVE_INTERPOLATOR_H_


#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>

//#include "cubic_spline.h"
#include "hashtable3_dh_class.h"
#include "cubic_spline.h"

namespace plant{

//typedef double Float;


/*
  0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16
  |           |           |           |           |
  |--.--|--.--|--.--|-^.--|--.--|--.--|--.--|--.--|
  |           |       !   |           |           |
  0           4       !   8           12          16      level 0 indices
  0     2     4     6 !   8     10    12    14    16      level 1 indices
  0  1  2  3  4  5  6 !7  8  9  10 11 12 13 14 15 16      level 2 indices = level max_depth indices = consecutive integers
                      !
                      !
  0     2     4       !   8  9  10 11 12          16      Indices stored in lookup table
  |-----|-----|-------^---|--|--|--|--|-----------|
                      !
                      +--- x = 6.7 (say). 
                          
Suppose SubdivisionLine looks like this ^^

                      Nearest index in level 0  = 4
                                    in level 1  = 6
                                    in level 2  = 6
                      
                      Since 6 is not stored in the LUT, 
                      but 4 is, x lies in a level 0 
                      interval starting at 4 

*/


#define __DEBUG_AI_ if(false)

inline int exp2i(unsigned int x){
	return 1 << x;
}


class SubdivisionSpline : public Spline{
	public:
	int npoints0 = 17;
	int max_depth = 16;
	double rel_tol = 1e-4;
	double abs_tol = 1e-4;
	
	//int npoints;
	int depth;
	int npoints_max;
	
	std::list <int> indices;
	
	HashTable<int,int> ht;
	
	double a, b, dx_min;
	std::list<double> xx, yy;
	std::list<bool> zz;

	public:
		
	SubdivisionSpline(){
	}


	SubdivisionSpline(int _n0, int _max_depth){
		npoints0 = _n0;
		depth = 0;
		max_depth = _max_depth;
		
		npoints_max = (npoints0-1)*exp2i(max_depth);
		dx_min = 1;
		a = 0; b = npoints_max;
		
		for (int i=0; i<npoints0-1; ++i){
			indices.push_back(i*npoints_max/(npoints0-1));
		}
		
	}
	
	
	// Subdivide interval [it, it+1] by inserting a point at the midpoint
	// Since the number of points is a power of 2, the midpoint will always be an integer
	std::list<int>::iterator split(std::list<int>::iterator& it){ // need the iterator by reference
		int a = *it;
		int b = *(++it);
		if (it == indices.end()){
			std::cout << "Already beyond the last interval" << std::endl;
			return indices.end();
		}
		if (b-a <= 1){
			std::cout << "Max refinement reached (a=" << a << ", b=" << b << ")\n";
			return indices.end();
		}
		else return indices.insert(it, (a+b)/2);	// integer division
	}
	
	
	void print() const{
		std::cout << "Subdivision Spline:\n";
		std::cout << "~~~~~~~~~~~~~~~~~~~\n";
//		for (auto i : indices) std::cout << i << " ";
//		std::cout << "| " << npoints_max << std::endl;
		std::cout << "depth = " << depth << ", maxdepth = " << max_depth << ", Nmax = " << npoints_max << std::endl;
		std::cout << a << " --- " << b << std::endl;
		std::cout << "Points (" << xx.size() << "):" << std::endl;
		auto xi = xx.begin();
		auto yi = yy.begin();
		auto zi = zz.begin();
		auto ii = indices.begin();
		auto mxi = m_x.begin(), myi = m_y.begin();
		for (; xi != xx.end(); ++xi, ++yi, ++zi, ++ii, ++mxi, ++myi) {
			std::cout << *ii << " (" << (b-a)*(*ii)/npoints_max << ") " << *xi << " " << *mxi << " " << *yi << " " << *myi << " | " << eval(*xi) << std::endl;
			assert(fabs((b-a)*(*ii)/npoints_max - *xi) < 1e-6);
		}
		std::cout << "~~~~~~~~~~~~~~~~~~~\n" << std::endl;
	}


	template <typename Function>  
	void construct(Function f, double _a, double _b){
		// auto t1 = std::chrono::steady_clock::now();
		__DEBUG_AI_ std::cout << "Constructing Spline: (" << _a << ", " << _b << ")\n";
		a = _a; b = _b;
		assert(a < b);
		
		npoints = npoints0;
		npoints_max = (npoints0-1)*exp2i(max_depth);
		depth = 0;
		dx_min = (b-a)/npoints_max;
		
		double dx = (b - a)/(npoints0 - 1);
		
		xx.clear(); yy.clear(); zz.clear(); indices.clear();
		for (int i = 0; i < npoints0; ++i) {
			const double xi = a + i*dx;
			xx.push_back(xi);
			yy.push_back(f(xi));
			zz.push_back(true);
			indices.push_back(i*npoints_max/(npoints0-1));
			
			__DEBUG_AI_ std::cout << xi << " " << f(xi) << " " << i*npoints_max/(npoints0-1) << std::endl;
			assert(!std::isnan(f(xi)));
		}
		__DEBUG_AI_ std::cout << "Size = " << xx.size() << " " << yy.size() << std::endl;
		
		set_points(xx, yy);

		bool b_error = true;
		while (b_error) {
			b_error = refine(f);
		}

		__DEBUG_AI_ std::cout << "Hashing intervals..." << std::endl;
		hashIntervals();
		__DEBUG_AI_ std::cout << "... DONE" << std::endl;
		
		// print();

		// auto t2 = std::chrono::steady_clock::now();
		// std::cout << "Interpolator Construction: (" << xx.size() << " points) [" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " usec]" << std::endl;
	}
	
	
	// Scale points
	// [a0, b0] ---> [_a, _b]
	template <typename Function>  
	void rescale(Function f, double _a, double _b){
		__DEBUG_AI_ std::cout << "Rescaling Spline: (" << _a << ", " << _b << ")\n";
		
		double a0 = xx.front();
		double b0 = xx.back();
		
		const double scale = (_b - _a) / (b0 - a0);
		for (auto itx = xx.begin(), ity = yy.begin(); itx != xx.end() && ity != yy.end(); ++itx, ++ity){
			double xnew = _a + (*itx - a0) * scale;
			double ynew = f(xnew);
			*itx = xnew;
			*ity = ynew;
		}
		set_points(xx,yy);
		
	}


	inline bool check_no_err(double y_true, double y_pred) const {
		const double err_abs = fabs(y_true - y_pred);
		const double err_rel = fabs(1 - y_pred / y_true);
		return err_abs < abs_tol || err_rel < rel_tol;
	}

	
	template <typename Function>
	bool refine(Function f){
		__DEBUG_AI_ std::cout << "refine... " << std::endl; 
		bool flag_refine = false;

		++depth;
		if (depth > max_depth){
			std::cout << "Maximum refinement reached around:" << std::endl;
			
		   	auto itx = std::next(xx.begin()), ity = std::next(yy.begin());
			for (auto it  = std::next(indices.begin());
					  it != std::prev(indices.end()); 
					  ++it, ++itx, ++ity){
				if ((*std::next(it) - *it) == 1 && (*it - *std::prev(it)) == 1){
					std::advance(it, -5);
					std::advance(itx, -5);
					std::advance(ity, -5);
					for (int i=0; i<10; ++i){
						std::cout << std::setw(10) << *it << " "
								  << std::setw(10) << *itx << " "
								  << std::setw(10) << *ity << "\n";
						++it; ++itx; ++ity;
					}
					std::cout << "---\n";
					std::advance(it, -5);
					std::advance(itx, -5);
					std::advance(ity, -5);
				}
			}
			print();

			//std::ofstream fout("interpolator_dump.txt");
			//for (int i=0; i<1000; ++i){
			//      double xm = 0.0998653;
			//      double xM = 0.0998729;
			//    double x = xm + (xM-xm)*i/1000.0;
			//    fout << x << " " << eval(x) << " " << f(x) << "\n";
			//}
			//fout.close();
			
			assert(false);
		}
		
		double dx = (b - a)/((npoints0 - 1)*exp2i(depth)); 
	
		auto xi = ++xx.begin(), yi = ++yy.begin();
		auto zi = ++zz.begin();
		auto ii = ++indices.begin();
		//std::cout << "depth = " << depth << std::endl;
		for (; xi != xx.end(); ++xi, ++yi, ++zi, ++ii) {
			if (*zi) {
				const double x_mid = *xi - dx;
				const double y_mid = f(x_mid);
				const double p_mid = eval(x_mid);
				__DEBUG_AI_ std::cout << "adding " << x_mid << ":\n"
									  << "   x:  " << std::setw(10) << *std::prev(xi)
									  << "       " << std::setw(10) << *xi - dx
									  << "       " << std::setw(10) << *(xi) << "\n"
									  << " f(x): " << std::setw(10) << f(*std::prev(xi))
									  << "       " << std::setw(10) << f(*xi - dx)
									  << "       " << std::setw(10) << f(*(xi))  << "\n"
									  << " i(x): " << std::setw(10) << eval(*std::prev(xi))
									  << "       " << std::setw(10) << eval(*xi - dx)
									  << "       " << std::setw(10) << eval(*(xi)) << "\n";

				// Always insert the new points (why waste the expensively computed true values!).
				xx.insert(xi, x_mid);
				yy.insert(yi, y_mid);
				split(--ii);

				// Check if this interval really needed refinement
				const bool flag_refine_mid = !check_no_err(y_mid, p_mid);

				// if so, mark both newly created intervals (on either side of the midpoint) for refinement
				*zi = flag_refine_mid;
				zz.insert(zi, flag_refine_mid);

				flag_refine = flag_refine || flag_refine_mid;
			}
		}

		// Recostruct spline to use new points added during refinement.
		set_points(xx, yy);

		return flag_refine;
	}


	// Below are functions to search for the interval in which x lies without
	//   binary searching the entire array

	void hashIntervals(){
		__DEBUG_AI_ std::cout << "assignment begins" << std::endl;
		ht = HashTable<int,int> (8093, 0.25);
		int ix  = 0;
		__DEBUG_AI_ std::cout << "hastable created" << std::endl;
		for (auto val : indices){
			__DEBUG_AI_ std::cout << "hashing " << val << std::endl;
			ht.hash_insert(val,ix);
			++ix;
		}
		//ht.hash_print();
	}


	mutable int nhteval_linear = 0;
	mutable int nhteval_bin = 0;
	mutable int n_lin = 0, n_bin = 0;
	

	// If x were to lie in an interval [a,b] in level 'level', find a. 
	int getIndexInLevel(double x, int level) const{
		int xi = x/dx_min;
		int dx = exp2i(max_depth-level); 
		return int(xi/dx)*dx;
	}


	// find the interval in which x lies 
	// by finding an interval at the highest depth which contains x 
	int find_interval(double x) const{
		++n_lin;
		//std::cout << "Find index linear:\n";
		int d = depth;
		while (d >= 0){
			//std::cout << "Depth " << d << std::endl;
			int idx = getIndexInLevel(x, d);
			auto h = ht.hash_find(idx);
			int it = h.id;
			nhteval_linear += h.attempts;

			if (it != -1){
				//std::cout << "~~ DONE\n";
				return  ht.ht[it].value; //distance(it,indices.begin());
			}
			--d;
		}
	}

	
	// This function uses the logic of finding first occurence from here:
	// https://stackoverflow.com/questions/38058740/binary-search-like-algorithm-to-find-change-in-value-in-sorted-array
	// It is modified slightly because the index runs in the opposite direction (depth-->0) 
	int find_interval_bin(double x) const{
		++n_bin;
		__DEBUG_AI_ std::cout << "Find index binary:\n";
		int result = 0;
		int start = depth;
		int end = 0;
		int it = -1;
		while(start >= end){
			int mid = (start + end)/2;
			//std::cout << "LargestIndex: " << start << " " << mid << " " << end << std::endl;
			__DEBUG_AI_ std::cout << "Depth " << mid << std::endl;
			int idx = getIndexInLevel(x, mid);
			auto h = ht.hash_find(idx);
			int it = h.id;
			nhteval_bin += h.attempts;
						
			if(it > -1){
				result = ht.ht[it].value;
				end = mid + 1;
			}
			else{
				start = mid - 1;
			}   
		}
		__DEBUG_AI_ std::cout << "~~ DONE\n";
		return result;
	}
	
	
	double eval_fast(double x) const{
		if (x < a) return extrapolate_left(x);
		else if (x >= b) return extrapolate_right(x);
		else {
			int interval_id = find_interval_bin(x);
			
			__DEBUG_AI_{			
				int id_linear = find_interval(x);
				std::cout << "compare liner and binary:" << interval_id << " " << id_linear << std::endl;
				assert(interval_id == id_linear);
			}
			
			__DEBUG_AI_ std::cout << "x-> " << x << " [" << m_x[interval_id] << ", " << m_x[interval_id+1] << "]" << std::endl;
			
			if (!(x >= m_x[interval_id]-1e-4 && x <= m_x[interval_id+1]+1e-4)){
				print();
				std::cout << "x-> " << x << " [" << m_x[interval_id] << ", " << m_x[interval_id+1] << "]" << std::endl;
			}
			
			// verify that correct interval was found
			assert(x >= m_x[interval_id]-1e-4 && x <= m_x[interval_id+1]+1e-4);
			
			return eval(x, interval_id);
		}
		
	}
	
};


typedef SubdivisionSpline AdaptiveInterpolator;

}

#endif

