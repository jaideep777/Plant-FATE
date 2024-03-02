#ifndef PSPM_CUBIC_SPLINE_H_
#define PSPM_CUBIC_SPLINE_H_

#include <vector>
#include <cassert>
#include <iomanip>
#include <limits>
typedef double Float;




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Thomas Algorithm to solve Tridiagonal system in O(N) time
//
//   n is the number of unknowns = matrix size
//   The system solved is as follows:
// 
//   b0 c0  0  0  0  0  0 	  x0     d0
//   a1 b1 c1  0  0  0  0     x1     d1
//    0 a2 b2 c2  0  0  0     x2     d2
//    0  0 a3 b3 c3  0  0  *  x3  =  d3
//    0  0  0 a4 b4 c4  0     x4     d4
//    0  0  0  0 a5 b5 c5     x5     d5
//    0  0  0  0  0 a6 b6     x6     d6
//
//   When storing the diagonal (b), upper (c), and lower (a) arrays, 
//   array index corresponds to row number in the matrix.  
//   Thus, a[0] and c[n-1] are never used.
//
//   Written by Keivan Moradi, 2014
//   See:
//   https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline void thomas_solve(double* a, double* b, double* c, double* d, int n) {
	n--; // since we start from x0 (not x1) JAI: Therefore last valid index is n for a,b, n-1 for c
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i]*c[i-1];
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i-- > 0;) {
		d[i] -= c[i]*d[i+1];		// JAI: note, i is in [0, n-1]
	}
}





// This is a custom implementation of std::lower_bound() for arrays
// The implementation can be found here:
//   http://www.cplusplus.com/reference/algorithm/lower_bound/
// This function returns an the index of the first element 
//   in the range [0,arrlen) which does not compare less than val.
//   i.e. index of first element >= val.
template <typename T>
int my_lower_bound(const T& val, const T * arr, const int arrlen){
	int it, first=0, last=arrlen;
	int count, step;
	count = arrlen;
	while (count>0){
		it = first; step=count/2; 
		it += step;
		if (arr[it] < val) {              
			first=++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}

// This is a custom implementation of std::upper_bound() for arrays
// The implementation can be found here:
//   http://www.cplusplus.com/reference/algorithm/upper_bound/
// This function returns an the index of the first element 
//   in the range [0,arrlen) which compares greater than val.
//   i.e. index of first element > val.
template <typename T>
int my_upper_bound(const T& val, T * arr, int arrlen){
	int it, first=0, last=arrlen;
	int count, step;
	count = arrlen;
	while (count>0){
		it = first; step=count/2; 
		it += step;
		if (!(val<arr[it])) { 
			first=++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}



// This spline interpolation follows the method from 
// https://kluge.in-chemnitz.de/opensource/spline/

// In interval [x{i} - x{i+1}], the function is 
// fi(x) = ai*(x-xi)^3 + bi*(x-xi)^2 + ci*(x-xi) + yi
// Applying constraints gives:
//  1.  ai = (b{i+1} - b{i})/3/hi     hi = x{i+1}-xi
//  2.  ci = (y{i+1}-y{i})/hi - 1/3*(2bi + b{i+1})hi
//  3.  Tridiagonal system for b: for i = [1,n-2]
//      h{i-1}b{i-1}/3 + 2(h{i-1}+h{i})bi/3 + hib{i+1}/3 = (y{i+1}-y{i})/hi - (y{i}-y{i-1})/h{i-1}
//  4.  Boundary conditions: (NEED TO BETTER UNDERSTAND)
//      

class Spline{
	public:
	int npoints;
	
	enum Type{LINEAR, CUBIC, CONSTRAINED_CUBIC};
	Type splineType = CUBIC;
	enum Extr{ZERO, CONSTANT, QUADRATIC, NA};
	Extr extrapolate = CONSTANT;

//	bool extrapolate_constant = true;
//	bool extrapolate_off = true;
	
	protected:
	std::vector <Float> m_a, m_b, m_c;
	std::vector <Float> m_x, m_y;		// Need random-access containers here
	
	public:
	template <typename Container>
	void set_points(const Container &x, const Container &y){	// construct should be able to take any containers, as the points will be copied to the Spline's own storage anyway
	   assert(x.size() == y.size());
	   m_x.assign(x.begin(), x.end());
	   m_y.assign(y.begin(), y.end());
	   
		npoints = x.size();
		m_a.resize(npoints);
		m_b.resize(npoints);
		m_c.resize(npoints);

	   // Check that x is in increasing order
	   for(int i=0; i<npoints-1; i++) {
		  assert(m_x[i]<m_x[i+1]);
	   }

	   if( splineType == CUBIC) { // cubic spline interpolation
			solve_coeffs_cubic();
	   } 
	   else if (splineType == LINEAR) { // linear interpolation
			solve_coeffs_linear();
	   }
	   else if (splineType == CONSTRAINED_CUBIC){
			solve_coeffs_constrained_cubic();	
	   }

	}


	inline void solve_coeffs_linear(){
		int n = npoints;  
		for(int i=0; i<n-1; i++) {
			 m_a[i]=0.0;
			 m_b[i]=0.0;
			 m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
		}
		m_a[n-1] = m_b[n-1] = 0.0;
		m_c[n-1] = m_c[n-2];	// fn-1 is same as fn-2 - same line extends right for extrapolation
	}


	inline void solve_coeffs_cubic(){
		int n = npoints;  
		  // setting up the matrix and right hand side of the equation system
		  // for the parameters b[]
		  //band_matrix A(n,1,1);
		  std::vector <Float> l(n,0), m(n,0), u(n,0);	// lower, middle, and upper diagonals of thr matrix. These hold coefficients of b{i-1}, b{i}, b{i+1} respectively
		  std::vector<double>  rhs(n);
		  for(int i=1; i<n-1; i++) {
			 l[i] =1.0/3.0*(m_x[i]-m_x[i-1]);
			 m[i] =2.0/3.0*(m_x[i+1]-m_x[i-1]);
			 u[i] =1.0/3.0*(m_x[i+1]-m_x[i]);
			 rhs[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]) - (m_y[i]-m_y[i-1])/(m_x[i]-m_x[i-1]);
		  }
		  // boundary conditions, zero curvature b[0]=b[n-1]=0
		  m[0]     = 2.0;
		  u[0]     = 0.0;
		  rhs[0]   = 0.0;	
		  m[n-1]   = 2.0;
		  l[n-1]   = 0.0;
		  rhs[n-1] = 0.0;	

		  // solve the equation system to obtain the parameters b[]
		  thomas_solve(l.data(),m.data(),u.data(),rhs.data(),n);
			m_b = rhs;

		  // calculate parameters a[] and c[] based on b[]
		  for(int i=0; i<n-1; i++) {
			 m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(m_x[i+1]-m_x[i]);
			 m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i])
					- 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(m_x[i+1]-m_x[i]);
		  }
		   // for the right boundary we define
		   // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
		   double h=m_x[n-1]-m_x[n-2];
		   // m_b[n-1] is determined by the boundary condition (turns out to be 0; = 0 for 0 curvature at xn-1)
		   m_a[n-1] = 0.0;
		   m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
	}


	inline void solve_coeffs_constrained_cubic(){
		int n = npoints;  

		std::vector <double> m(n);
		for (int i=1; i<n-1; ++i){
			double fplus  = (m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
			double fminus = (m_y[i]-m_y[i-1])/(m_x[i]-m_x[i-1]);
			if (fplus*fminus < 0) m[i] = 0;
			else m[i] = 2/(1/fplus + 1/fminus);  // harmonic mean will tend towards smaller derivative
			//else m[i] = (fplus+fminus)/2; 
		}
		// Boundary conditions: zero curvature at x0 and xN-1 
		m[0] = 3*(m_y[1]-m_y[0])/(m_x[1]-m_x[0])/2 - m[1]/2;
		m[n-1] = 3*(m_y[n-1]-m_y[n-2])/(m_x[n-1]-m_x[n-2])/2 - m[n-2]/2;

		for (int i=0; i<n-1; ++i){
			m_c[i] = m[i];
			double hi = m_x[i+1]-m_x[i];
			m_b[i] = ( 3*(m_y[i+1]-m_y[i] - m[i]*hi) - hi*(m[i+1]-m[i]))/hi/hi;
			m_a[i] = (-2*(m_y[i+1]-m_y[i] - m[i]*hi) + hi*(m[i+1]-m[i]))/hi/hi/hi;
		}
		m_b[n-1] = 0;
		m_c[n-1] = m[n-1];
	}


	template <typename Container>
	void constructAndReset(const Container &x, const Container &y){	// construct should be able to take any containers, as the points will be copied to the Spline's own storage anyway
		set_points(x,y);
	}

	inline Float extrapolate_left(double x) const{
		  double h=x-m_x[0];
		  // extrapolation to the left
		  if (extrapolate == ZERO) return 0;	// zero
		  else if (extrapolate == NA) return std::numeric_limits<double>::quiet_NaN();	// nan
		  else if (extrapolate == CONSTANT) return m_y[0];	// constant
		  else return ((m_b[0])*h + m_c[0])*h + m_y[0];  // quadratic
	}
	
	inline Float extrapolate_right(double x) const{
			size_t n=m_x.size();
			double h=x-m_x[n-1];
		  // extrapolation to the right
			if (extrapolate == ZERO) return 0;	// zero
			else if (extrapolate == CONSTANT) return m_y[n-1]; // constant
			else return ((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];  // quadratic
	
	}
	
	inline Float eval(Float x) const {
		size_t n=m_x.size();

		if(x<m_x[0]) return extrapolate_left(x);
		else if(x>=m_x[n-1]) return extrapolate_right(x);
		else {
			// find the closest point m_x[idx] < x
//			std::vector<double>::const_iterator it;
//			it=std::lower_bound(m_x.begin(),m_x.end(),x);
//			int idx=std::max( int(it-m_x.begin())-1, 0);
			int it = my_lower_bound(x, m_x.data(), m_x.size());
			int idx=std::max(it-1,0); //std::max(it-1, 0);
//			int it_ub = my_upper_bound(x, m_x.data(), m_x.size());
//			int idx_ub=it_ub-1; //std::max(it_ub-1, 0);
//			std::cout << "x: " << x << " " << idx << " " << idx_ub << endl;
			
			return eval(x,idx);			
		}
	}
	
	inline Float eval(Float x, int idx) const{
		double h=x-m_x[idx];
		return  ((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
	}

	inline void print(){
		std::cout << "Spline: \n";
		for (int i=0; i<npoints; ++i){
			std::cout << std::right << std::setw(10) << std::setfill(' ') << m_x[i] << " " 
					  << std::right << std::setw(10) << std::setfill(' ') << m_y[i] << " "
					  << std::right << std::setw(10) << std::setfill(' ') << m_a[i] << " "
					  << std::right << std::setw(10) << std::setfill(' ') << m_b[i] << " "
					  << std::right << std::setw(10) << std::setfill(' ') << m_c[i] << " " 
					  << "\n";
		}
		std:: cout << std::endl;
	}
};



// When storing the diagonal, upper, and lower arrays, array
// index corresponds to row number. Thus, a[0] and c[N-1] are 
// never used
 
class BandMatrix{
	public:
	int N;
	
	std::vector <Float> upper;	// element N-1 is absent
	std::vector <Float> diag;
	std::vector <Float> lower;	// element 0 is absent

	inline BandMatrix(int size){
		N = size;
		upper.resize(N);
		lower.resize(N);
		diag.resize(N);
	}

	inline BandMatrix(std::vector<Float> _lower, std::vector<Float> _diag, std::vector<Float> _upper){
		N = _upper.size();
		upper = _upper;
		lower = _lower;
		diag  = _diag;
	}

	inline void print(){
		for (int i=0; i<N; ++i){
			for (int z=0; z<i-1;++z) std::cout << 0 << "\t";
			if (i > 0) std::cout << lower[i] << "\t";
			std::cout << diag[i] << "\t";
			if (i < N-1) std::cout << upper[i] << "\t";
			for (int z=i+2; z<N; ++z) std::cout << 0 << "\t";
			std::cout << "\n";
		}
	}
	
	inline void solve(std::vector<Float>&y){
		thomas_solve(lower.data(), diag.data(), upper.data(), y.data(), N);
	}
};



#endif

