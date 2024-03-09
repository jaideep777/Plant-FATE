#include <iomanip>
#include <fstream>
#include <random>

#include "life_history.h"
using namespace std;

int main(){

	cout << setprecision(12);
	
	pfate::LifeHistoryOptimizer lho("tests/params/p_test_v2.ini");
	// lho.C.init_co2(414);

	// lho.ts.set_units("days since 0000-01-00 0:0:0");
	lho.ts.set_units("years since 0000-01-00 0:0:0");
	lho.par0.set_tscale(lho.ts.get_tscale());

	lho.init();
	lho.C.Climate::print(0);

	cout << "Starting rates:\n";
	lho.P.calc_demographic_rates(lho.C, 0);
	cout << "  npp   = " << lho.P.res.npp   *365 / lho.par0.days_per_tunit << '\n';
	cout << "  rleaf = " << lho.P.res.rleaf *365 / lho.par0.days_per_tunit << '\n';
	cout << "  rstem = " << lho.P.res.rstem *365 / lho.par0.days_per_tunit << '\n';
	cout << "  rroot = " << lho.P.res.rroot *365 / lho.par0.days_per_tunit << '\n';
	cout << "  tleaf = " << lho.P.res.tleaf *365 / lho.par0.days_per_tunit << '\n';
	cout << "  troot = " << lho.P.res.troot *365 / lho.par0.days_per_tunit << '\n';
	cout << "  mort  = " << lho.P.rates.dmort_dt *365 / lho.par0.days_per_tunit << '\n';
	cout << "  fec   = " << lho.P.rates.dseeds_dt *365 / lho.par0.days_per_tunit << '\n';

	if (fabs(lho.P.res.npp   *365 / lho.par0.days_per_tunit / 0.00649373522273 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.res.rleaf *365 / lho.par0.days_per_tunit / 0.00105810720763 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.res.rstem *365 / lho.par0.days_per_tunit / 0.00223950139619 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.res.rroot *365 / lho.par0.days_per_tunit / 0.00063537663421 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.res.tleaf *365 / lho.par0.days_per_tunit / 0.000356239104633 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.res.troot *365 / lho.par0.days_per_tunit / 0.000962086104036 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.rates.dmort_dt *365 / lho.par0.days_per_tunit / 0.0658361385886 - 1) > 1e-6) return 1; 
	if (fabs(lho.P.rates.dseeds_dt *365 / lho.par0.days_per_tunit / 0.000402691681301 - 1) > 1e-6) return 1; 
	// if (fabs(lho.P.rates.dmort_dt / 2.65329323718 - 1) > 1e0.0658361385886 return 1; 

	double total_prod = lho.P.get_biomass();
	cout << "Starting biomass = " << total_prod << "\n";
	cout << "Mortality until seedling stage = " << lho.P.state.mortality << "\n";

	double dpt = 365.2425/lho.par0.days_per_tunit;

	ofstream fout("assim1.txt");
	double dt = 0.1*dpt;
	lho.printHeader(fout);
	for (double t=2000*dpt; t<=2500*dpt; t=t+dt){
		lho.grow_for_dt(t, dt);
		lho.printState(t+dt, fout);
		// lho.C.Climate::print(t);
	}
	fout.close();
	// fswp.close();

	cout << "At last t: " << "\n" 
		 << "  Total biomass    = " << lho.P.get_biomass() << "\n"
		 << "  Total litter     = " << lho.litter_pool << "\n"
		 << "  Total reproduc   = " << lho.rep << "\n"
		 << "  Total bio+lit+rep = " << lho.P.get_biomass() + lho.litter_pool + lho.rep << "\n"
		 << "  Total production = " << lho.prod << "\n";
	

	double rel_error = abs((lho.P.get_biomass()+lho.litter_pool+lho.rep)/lho.prod - 1);
	cout << "Relative error in biomass accounting = " << rel_error << endl;
	if (rel_error > 2e-5) return 1;
	
	double fitness1 = lho.seeds;
	cout << "Fitness = " << fitness1 << '\n';
	// if (fabs(fitness1 - 0.0999331152132) > 1e-6) return 1;  // expected value updated after implementing inst phydro and setting kphio to 0.045
	// if (fabs(fitness1 - 0.0998537505871) > 1e-6) return 1;  // expected value updated after using actual p88/p50 ratio instead of 3.01
	// if (fabs(fitness1 - 0.0192874440426) > 1e-6) return 1;  // expected value updated after using daylength in phydro
	if (fabs(fitness1 - 0.0192813597174) > 1e-6) return 1;  // expected value updated after bugfix in unit conversion in assimilation.tpp


	// FIXME. This reinit of LHO does not work. assim rate is NAN.
	// lho.init();
	// lho.C.Climate::print(0);
	// total_prod = lho.P.get_biomass();
	// cout << "Starting biomass = " << total_prod << "\n";
	// cout << "Mortality until seedling stage = " << lho.P.state.mortality << "\n";
	// double fitness = lho.calcFitness();
	// cout << "After calcFitness(): " << "\n" 
	// 	 << setprecision(12) 
	// 	 << "  Total biomass    = " << lho.P.get_biomass() << "\n"
	// 	 << "  Total litter     = " << lho.litter_pool << "\n"
	// 	 << "  Total reproduc   = " << lho.rep << "\n"
	// 	 << "  Total bio+lit+rep = " << lho.P.get_biomass() + lho.litter_pool + lho.rep << "\n"
	// 	 << "  Total production = " << lho.prod << "\n";
	// cout << "Fitness = " << fitness << endl;

	// if (fabs(fitness - 0.414567339728) > 1e-6) return 1;
	// if (fabs(fitness - 0.427527753304) > 1e-6) return 1; // expected value after upgrade to latest version of phydro @6fc30d6
	// if (fabs(fitness - 0.406735962511) > 1e-6) return 1; // expected value updated after finding minor bug in gpp calc... when applying midday --> day mean conversion, it was directly applied to gpp, whereas it should only be applied to a, since vcmax is not scaled 
	// if (fabs(fitness - 0.253198528939) > 1e-6) return 1;    // expected value updated after implementing inst phydro and setting kphio to 0.045


	return 0;
}

//# R script to test:
//dat = read.delim("/home/jaideep/codes/plant_fate_ppa/assim.txt")
//dat$heartwood_fraction = 1-dat$sapwood_fraction

//par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

//plot(dat$height~dat$diameter)
//plot(dat$crown_area~I(dat$height*dat$diameter))
//plot(dat$crown_area~dat$height)
//plot(dat$crown_area~dat$diameter)
//plot(dat$leaf_area~dat$crown_area)
//plot(dat$sapwood_fraction~dat$height)
//# plot(sqrt(4*dat$crown_area/pi)~dat$height)

//plot(dat$height~dat$i)
//plot(dat$diameter~dat$i)
//plot(dat$total_mass~dat$i)
//points(dat$total_prod~dat$i, type="l", col="red")

//matplot(y=cbind(dat$assim_gross,
//                dat$assim_net), 
//        x=dat$i, col=c("green3", "green4"), log="", pch=20)

//matplot(y=cbind(dat$rr,
//                dat$rs,
//                dat$rl), 
//        x=dat$i, col=c("yellow2", "yellow3", "yellow4"), log="", pch=20)

//matplot(y=cbind(dat$tr,
//                dat$tl), 
//        x=dat$i, col=c("orange3", "orange4"), log="", pch=20)


//D = dat$diameter
//H = dat$height
//-log(1-H/20)*20/D
