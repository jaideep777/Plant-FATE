
/************************************************************************
 *		z					C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;				a range [ax,bx]
 *	double f(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *			 the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************/
 

#include <cmath>
#include <limits>

namespace pn{

struct BrentRes{
	double root;
	double nfnct;
	BrentRes(double r, double n){root=r; nfnct=n;}
};

template <typename Function>
BrentRes zero(		// An estimate to the root	
	double ax,		// Left border | of the range	
	double bx,		// Right border| the root is seeked
	Function f,		// Function under investigation	
	double tol){	// Acceptable tolerance		

	double a,b,c;			// Abscissae, descr. see above	
	double fa;				// f(a)				
	double fb;				// f(b)				
	double fc;				// f(c)				

	int ncalls = 0;
	
	a = ax;	b = bx;	fa = f(a);	fb = f(b);
	c = a;	 fc = fa;
	ncalls += 2;

	for(;;){		// Main iteration loop	
	
		double prev_step = b-a;	// Distance from the last but one
							// to the last approximation	
		double tol_act;			// Actual tolerance		
		double p;						// Interpolation step is calcu- 
		double q;						// lated in the form p/q; divi- 
								// sion operations is delayed	 
 							// until the last moment	
		double new_step;				// Step at this iteration			 
	 
		if( fabs(fc) < fabs(fb) ){
			a = b;	b = c;	c = a;					// Swap data for b to be the best approximation	 		
			fa=fb;	fb=fc;	fc=fa;
		}
		tol_act = 2*std::numeric_limits<double>::epsilon()*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if( fabs(new_step) <= tol_act || fb == (double)0 ) return BrentRes(b, ncalls);	// Acceptable approx. is found	

					// Decide if the interpolation can be tried	
		if( fabs(prev_step) >= tol_act		// If prev_step was large enough
			&& fabs(fa) > fabs(fb) )	// and was in true direction,	
		{					// Interpolatiom may be tried	
			register double t1,cb,t2;
			cb = c-b;
			if( a==c )			// If we have only two distinct	
			{					// points linear interpolation 	
				t1 = fb/fa;		// can only be applied		
				p = cb*t1;
				q = 1.0 - t1;
		 	}
			else				// Quadric inverse interpolation
			{
				q = fa/fc;	t1 = fb/fc;	t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if( p>(double)0 )	// p was calculated with the op-
				q = -q;			// posite sign; make p positive	
			else				// and assign possible minus to	
				p = -p;			// q				

			if(		p < (0.75*cb*q-fabs(tol_act*q)/2)	// If b+p/q falls in [b,c]
				&& p < fabs(prev_step*q/2) )			// and isn't too large	
				new_step = p/q;			// it is accepted	
							// If p/q is too large then the	
							// bissection procedure can 	
							// reduce [b,c] range to more	
							// extent			
		}

		if( fabs(new_step) < tol_act )	// Adjust the step to be not less than tolerance 
			if( new_step > (double)0 ) new_step = tol_act;
			else new_step = -tol_act;

		a = b;	fa = fb;			// Save the previous approx.	
		b += new_step;	fb = f(b);	// Do step to a new approxim.	
		ncalls += 1;
		
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ){								 			
			c = a;	fc = fa;		// Adjust c for it to have a sign opposite to that of b				
		}
	}

}


}	// pn 
