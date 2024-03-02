#include <solver.h>
#include "test_model_2_ms.h"
using namespace std;

int main(){
	TestModel M;
	Environment E;

	Solver S(SOLVER_EBT);
	S.control.cm_use_log_densities = true;
	S.setEnvironment(&E);

	Species<TestModel> spp(M);

	//Solver<TestModel> S({0,1,2,3}, SOLVER_CM);
	// S.addSpecies({0,1,2,3}, &spp, 4, 2);
	S.addSpecies(vector<vector<double>>{{0,1,2,3}}, (Species_Base*)&spp, 4, 2);

	S.print();

	//     x
	//     ^ 
	//     |
	//   3 |---------------* (8)
	//     |  A = 6
	//   2 |-------*         (4)
	//     |  A = 3
	//   1 |---*             (2)
	//     |  A = 1.5 
	//   0 |-*               (1)
	//     +-----------------------> u

	S.state = {3,(8), 2,(4), 1,(2), 0,(0), 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
	S.copyStateToCohorts(S.state.begin());
	S.print();

	auto w = [](int i, double t){return 1;};

	cout << "State integral = " << S.state_integral(w, 0, 0) << endl; //<< " " << S.integrate_wudx_above(w, 0, 0, 0) << endl;
	
	if(abs(S.state_integral(w, 0, 0) - 14) > 1e-5) return 1;
	// if(abs(S.state_integral(w, 0, 0) - S.integrate_wudx_above(w, 0, 0, 0)) > 1e-5) return 1;
	// if(abs(S.integrate_wudx_above(w, 0, 1.5, 0) - 7.75) < 2e-5) return 1;


	// //     x
	// //     ^ 
	// //     |
	// //   4 |---------------* (8)
	// //     |  A = 6
	// //   3 |-------*         (4)
	// //     |  A = 3
	// //   2 |---*             (2)
	// //     |  A = 1.5 
	// //   1 |-*               (1)
	// //	  |  A = 1.5   
	// //   0 +----*------------(2)---> u
	// //
	// S.state = {4,log(8), 3,log(4), 2,log(2), 1,log(1), 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
	// S.copyStateToCohorts(S.state.begin());
	// S.control.integral_interpolate = false;
	// //cout << S.state_integral(w, 0, 0) << " " << S.integrate_wudx_above(w, 0, 1.1, 0) << endl;
	// //if(abs(S.state_integral(w, 0, 0) - 12) > 1e-5) return 1;
	// //if(abs(S.integrate_wudx_above(w, 0, 1.1, 0) - 10.5) > 1e-5) return 1;

	// vector<double> xlow = {0,0.5, 1,1.5, 2};
	// vector<double> expected = {12, 12, 10.5, 10.5, 9};
	// for (int i=0; i<xlow.size(); ++i){
	// 	double I = S.integrate_wudx_above(w,0,xlow[i],0);
	// 	cout << "I(" << xlow[i] << ") = " << I << "\n";
	// 	if(abs(I - expected[i]) > 1e-5)	return 1;
	// }
	// cout << "Ix(0) = " << S.state_integral(w, 0, 0) << " " << endl;
	// if(fabs(S.state_integral(w, 0, 0) - 12) > 1e-5) return 1; 
	

	// //     x
	// //     ^ 
	// //     |
	// //   4 |---------------* (8)
	// //     |  A = 6
	// //   3 |-------*         (4)
	// //     |  A = 3
	// //   2 |---*             (2)
	// // 1.1 |  *   A = 1.5    (1.1)
	// //   1 |-*               (1)
	// //	   |  A = 1.5   
	// //   0 +----*------------(2)---> u
	// //
	// S.control.integral_interpolate = true;
	// xlow = {0,0.5, 1,1.5, 2};
	// expected = {12, 11.125, 10.5, 9.875, 9};
	// for (int i=0; i<xlow.size(); ++i){
	// 	double I = S.integrate_wudx_above(w,0,xlow[i],0);
	// 	cout << "I(" << xlow[i] << ") = " << I << "\n";
	// 	if(abs(I - expected[i]) > 1e-5)	return 1;
	// }
	// cout << "Ix(0) = " << S.state_integral(w, 0, 0) << " " << endl;
	// if(fabs(S.state_integral(w, 0, 0) - 12) > 1e-5) return 1; 

	// //S.state = {4,log(8), 3,log(4), 2,log(2), 1,log(1), 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
	// //S.copyStateToCohorts(S.state.begin());
	// //S.control.integral_interpolate = true;
	// //cout << S.state_integral(w, 0, 0) << " " << S.integrate_wudx_above(w, 0, 1.1, 0) << endl;
	// //if(abs(S.state_integral(w, 0, 0) - 12) > 1e-5) return 1;
	// //if(abs(S.integrate_wudx_above(w, 0, 1.1, 0) - 10.395) > 1e-5) return 1;


	return 0;

}
