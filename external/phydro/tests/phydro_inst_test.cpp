#include "phydro.h"
using namespace std;
using namespace phydro;

double err(double x, double ref){
	return std::min(abs(x-ref), abs(x/ref-1));
}

int check(double x, double ref, double err=1e-4){
	//cout << "err: " << abs(x-ref) << " " << abs(x/ref-1)<< "\n";
	//cout << "comp: " << x << " " << ref << "|" << abs(x-ref) << "\n";
	if (abs(x-ref) < err || abs(x/ref-1) < err) return 0;
	else return 1;
}


vector<double> seq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(start + double(i)*(end-start)/(length-1));
	return x;
}

int main(){
	
	double tc = 20;
	double vpd = 810.6;
	double co2 = 400;
	double ppfd = 1200;
	double kphio = 0.087;
	double fapar = 0.99;
	double rdark = 0.02;
	double psi_soil = -0.4137931;
	double elv = 0;
	double pa = phydro::calc_patm(elv);
	
	ParCost        par_cost(0.118514, 1.227068);
	ParPlant       par_plant(7.457324e-17, -1.039539, 1);
	ParPhotosynth  par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark, tc, tc);
	ParEnv         par_env(tc, pa, vpd, ppfd/2);

	double jmax = 117.0184518;
	double vcmax = 55.6279401;

	PHydro_Profit_Inst P(vcmax, jmax, psi_soil, par_cost, par_photosynth, par_plant, par_env);
	for (int i=0; i<50; ++i){
		double dpsi = 2.0/49.0*i;
		VectorXd x(1);
		x << dpsi;
		double grad = calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost);
		cout << dpsi << " " << P.value(x) << " " << grad << "\n";
	}

	double dpsi_opt_num = optimize_shortterm_multi(vcmax, jmax, psi_soil, par_cost, par_photosynth, par_plant, par_env);
	cout << "Opt = " << dpsi_opt_num << "\n";

	auto dpsi_opt = pn::zero(0, 20, [&](double dpsi){return calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost);}, 1e-6);
	cout << "Opt analytical = " << dpsi_opt.root << "\n";

	for (auto psi_soil : seq(-6, 0, 20)){

		phydro::PHydroResult res;
		res = phydro::phydro_instantaneous_analytical(vcmax, jmax, tc, tc, ppfd, ppfd/2, vpd, co2, pa, fapar, kphio, psi_soil, rdark, 3.0, par_plant, par_cost);
			
		cout << setw(10) <<  psi_soil  << "\t"; cout.flush();
		cout << setw(10) <<  res.jmax  << "\t";
		cout << setw(10) <<  res.dpsi  << "\t";
		cout << setw(10) <<  res.gs    << "\t";
		cout << setw(10) <<  res.a     << "\t";
		cout << setw(10) <<  res.ci    << "\t";
		cout << setw(10) <<  res.chi   << "\t";
		cout << setw(10) <<  res.vcmax << "\n"; cout.flush();
	}

	if (fabs(dpsi_opt_num - dpsi_opt.root) < 1e-5) return 0;
	else return 1;

	return 0;
}

