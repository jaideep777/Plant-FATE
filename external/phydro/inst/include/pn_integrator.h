// -*-c++-*-
#ifndef PLANT_PLANT_PN_INTEGRATOR_H_
#define PLANT_PLANT_PN_INTEGRATOR_H_

#include <cassert>
#include <chrono>
#include <cmath>

namespace pn {

class Integrator{
	private:
	
	size_t limit = 1000;
	double rel_tol = 1e-4;
	double abs_tol = 1e-4;
	
	public:
	double gauss = 0;
	double kronrod = 0;
	double err = 0;
	
	public:
	Integrator(){
//		gslws = gsl_integration_workspace_alloc(limit);
	}
	
	~Integrator(){
//		gsl_integration_workspace_free(gslws);
	}
	
	template <typename Func>  
	double integrate(Func &f, double a, double b){
		return qng21(f,a,b);
	}
	
	
	template <typename Func>  
	double qng21(Func &f, double a, double b){
		// QNG integrator mine
		// * QK21:
		// Abscissae of the 21-point kronrod rule:
		// xgk[1], xgk[3], ... abscissae of the 10-point gauss rule.
		// xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule
		const double xgk[11] = {
		  0.995657163025808080735527280689003, //  K
		  0.973906528517171720077964012084452, // G
		  0.930157491355708226001207180059508, //  K
		  0.865063366688984510732096688423493, // G
		  0.780817726586416897063717578345042, //  K
		  0.679409568299024406234327365114874, // G
		  0.562757134668604683339000099272694, //  K
		  0.433395394129247190799265943165784, // G
		  0.294392862701460198131126603103866, //  K
		  0.148874338981631210884826001129720, // G
		  0.000000000000000000000000000000000  //  K
		};

		// Weights of the 10-point gauss rule:
		const double wg[5] = {
		  0.066671344308688137593568809893332, // G
		  0.149451349150580593145776339657697, // G
		  0.219086362515982043995534934228163, // G
		  0.269266719309996355091226921569469, // G 
		  0.295524224714752870173892994651338  // G
		};

		// Weights of the 21-point kronrod rule:
		const double wgk[11] = {
		  0.011694638867371874278064396062192, //  K
		  0.032558162307964727478818972459390, // G
		  0.054755896574351996031381300244580, //  K
		  0.075039674810919952767043140916190, // G
		  0.093125454583697605535065465083366, //  K
		  0.109387158802297641899210590325805, // G
		  0.123491976262065851077958109831074, //  K
		  0.134709217311473325928054001771707, // G
		  0.142775938577060080797094273138717, //  K
		  0.147739104901338491374841515972068, // G
		  0.149445554002916905664936468389821  //  K
		};

		double G = 0, K = 0;
		auto t = [a, b](double x){
/*			cout << "x: " << x << " " << x_min + (x_max - x_min)/2*(x - -1) << endl;*/
			return a + (b - a)/2*(x + 1);
		};
		for (int i=0; i<5; ++i){
			int igk = 2*i+1;	// odd indices are common (Gauss and Kronrod)
			double fplus  = f(t(xgk[igk]));
			double fminus = f(t(-xgk[igk]));
			G += (fplus + fminus)*wg[i];
			K += (fplus + fminus)*wgk[igk];
		}
		for (int i=0; i<5; ++i){
			int ik = 2*i;	// Even indices are Kronrod
			K += (f(t(xgk[ik])) + f(t(-xgk[ik])))*wgk[ik];
		}
		K += f(t(xgk[10]))*wgk[10];
		
		G *= (b-a)/2;
		K *= (b-a)/2;

		gauss = G;
		kronrod = K;
		err = fabs(G-K);
		
		return K;
	}
	
};


template <typename ContainerX, typename ContainerY>
double integrate_trapezium(const ContainerX& x, const ContainerY& y) {
  assert(y.size() == x.size());
  assert(x.size() >= 2);
  
  auto x0 = x.begin(), x1 = ++x.begin();
  auto y0 = y.begin(), y1 = ++y.begin();
  double tot = 0.0;
  while (x1 != x.end()) {
    tot += (*x1++ - *x0++) * (*y1++ + *y0++);
  }
  return tot * 0.5;
}


}


#endif

