#include "ode_solver.h"
#include <fstream>
using namespace std;

int nsys = 0;

void sys(double x, vector<double>::iterator _y, vector<double>::iterator _dydx, void* _data){ 
  double* y = &*_y;
  double* dydx = &*_dydx;
  // Exact solution is y = log(exp(y0)+20t)
  for (int i=0; i<nsys; ++i) dydx[i] = 20*exp(-y[i]/1);
}

void after_step(double x, vector<double>::iterator _y){
    double* y = &*_y;
	cout << "After step: t = " << x << "\n";
}


int main(){
  for (auto str : {"lsoda", "rk45ck"}){
	  double t_start=0, t_stop=10; // elapsed time in dimensionless units

	  nsys = 0;
	  vector<double> y;

	  OdeSolver stepper(str, t_start, 1e-8, 1e-8); // RK class for adaptive step
	  
	  //cout.precision(15);
	  
	  double ti = t_start;
	  
	  ofstream fout((string("sc_sys_")+str+".txt").c_str());
	  for (double t=t_start; t <= t_stop; t=t+0.1) {
		//cout << t; // << " " << int(t);
		fout << ti << "\t";
		for (auto yy: y) fout << yy << "\t"; //[0]<<" "<<y[1]<<" | "<<y[0]-xmax*sin(ti)<<" | "<<y[1] - xmax*cos(ti)<< " | "<< Energy(y)<<endl;
		fout << "\n";
		
		stepper.step_to(t, ti, y,  sys, after_step);
		
		if (abs(t-int(t*1.00000001)) < 1e-6){ // insert cohort at integer t
			//cout << " insert.";
			nsys += 1;
			y.resize(y.size()+1, 0.1);
		}    
		//cout << "\n";
		

	  }
	  fout.close();
	 
		cout << "Number of fn evaluations " << str << " = " << stepper.get_fn_evals() << endl;  
	  //if (rk.size() != 2) return 1;
  }

  return 0;
}

//# plot with R:
//dat = read.delim("/home/jaideep/codes/pspm_ode_solver/sc_sys_lsoda.txt", header = F, sep = "\t", col.names = paste("V", 0:10))
//matplot(x=dat$V.0, y=dat[,-1], type="l", lty=1, col=rainbow(12))
//dat = read.delim("/home/jaideep/codes/pspm_ode_solver/sc_sys_rk45ck.txt", header = F, sep = "\t", col.names = paste("V", 0:10))
//matplot(x=dat$V.0, y=dat[,-1], type="l", lty=1, col=rainbow(12))
//points(x=seq(0,10,1), y=log(exp(0.1)+20*seq(0,10,1)))

