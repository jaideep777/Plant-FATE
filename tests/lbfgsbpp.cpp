#include <iostream>
#include <Eigen/Core>
using namespace std;

#include <LBFGSB.h>

using Eigen::VectorXd;
using namespace LBFGSpp;

class Rosenbrock{
	private:
    int n;
	double scale;

	public:
    Rosenbrock(int n_, double sc) : n(n_), scale(sc) {}

	double value(const VectorXd &x) {
		//const double t1 = (1 - x[0]);
		//const double t2 = (x[1] - x[0] * x[0]);
		//return   t1 * t1 + scale * t2 * t2;
		return  (x[0]-2)*(x[0]-2) + (x[1]-3)*(x[1]-3);
	}

    double operator()(const VectorXd& x, VectorXd& grad){
		double f = value(x);
		
		for (int i=0; i<n; ++i){
			VectorXd dx = VectorXd::Zero(n);
			double delta = 2.2204e-6;
			dx[i] = delta;
			double fplus = value(x+dx);
			double fminus = value(x-dx);
			grad[i] = (fplus-fminus)/delta/2;
		}
	
		return f;
    }
};


int main(){
    const int n = 2;
    // Set up parameters
    LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);
    Rosenbrock fun(n, 100);

	// bounds
	VectorXd lb(n), ub(n);
	lb << 0, 0;
	ub << 1, 5;


    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

	if (abs(x[0]-1) < 1e-6 && abs(x[1]-3) < 1e-6) return 0;
	else return 1;
}




