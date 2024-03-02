#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "daphnia.h"

int main(){

	Species<Daphnia> spp;
	Environment E;

	Solver S(SOLVER_FMU);
	S.control.ode_ifmu_stepsize = 0.1;
	
	S.setEnvironment(&E);

	S.addSpecies({300}, {0}, {1}, {false}, &spp, 0, -1);
	S.addSystemVariables({E.K});  // this can be done either before or after addSpecies()

	S.control.ebt_ucut = 1e-20;

	S.initialize();
	//S.print();
	
	ofstream fout("fmu_Daphnia.txt");

	for (double t=0.05; t <= 100; t=t+1) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
		cout << S.current_time << " " << S.state[0] << "\n";
		for (auto y : S.state) fout << y << "\t";
		
		fout << endl;
	}

	fout.close();

	// Expected 44.3530812 
	cout << S.newborns_out(100)[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

//# R script to verify code:
//dat = read.delim("/home/jaideep/codes/libpspm/demo/Daphina_model/fmu_Daphnia.txt", header=F)
//a = 0.75
//mu = 0.1
//r = 0.5
//K = 3
//xstar = (mu*(1+mu)*(2+mu)/2/a)^(1/3)
//sstar = xstar/(1-xstar)
//N = 300
//x=seq(0,1,length.out = N)
//plot(x=x, y=exp(-8*x^3), type="l")

//dat = dat[,-ncol(dat)]
//matplot(x=seq(0,1,length.out = N), y = t(dat[,-c(1,2,3)]), col=rainbow(100, start=0, end=0.9, alpha=0.2), type="l", lty=1, log="y", ylim=c(0.1, 1e4))
//abline(v=xstar, col="black")
//xeq = seq(0,1,length.out = 30)
//points(x=xeq, y= a*r*sstar*(1-sstar/K)*(xstar-xeq)^(mu-1)/xstar^mu)
//plot(dat$V2~dat$V1, type="l")




