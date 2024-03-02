#include <iostream>
#include <Eigen/Core>
using namespace std;

#include <LBFGS.h>

using Eigen::VectorXd;
using namespace LBFGSpp;

class Rosenbrock{
	private:
    int n;
	
	public:
    Rosenbrock(int n_) : n(n_) {}

	double value(const VectorXd &x) {
		const double t1 = (1 - x[0]);
		const double t2 = (x[1] - x[0] * x[0]);
		return   t1 * t1 + 100 * t2 * t2;
	}

    double operator()(const VectorXd& x, VectorXd& grad)
    {
		//double f = value(x);
		
		//for (int i=0; i<n; ++i){
		//    VectorXd dx = VectorXd::Zero(n);
		//    dx[i] = 1e-6;
		//    double fplus = value(x+dx);
		//    grad[i] = (fplus-f)/1e-6;
		//}
	
		//return f;

	    double fx = 0.0;
		for(int i = 0; i < n; i += 2)
		{
			double t1 = 1.0 - x[i];
			double t2 = 10 * (x[i + 1] - x[i] * x[i]);
			grad[i + 1] = 20 * t2;
			grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
    }
};


int main(){
    const int n = 2;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Rosenbrock fun(n);

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
