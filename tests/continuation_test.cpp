#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "plantfate.h"

using namespace std;

int main(int argc, char ** argv){

	string par_file = "tests/params/p_cont_test.ini";
	if (argc == 2){
		par_file = argv[1];
	}

	// reference run
	{
		Simulator sim(par_file);
		sim.continuePrevious = false;
		sim.expt_dir = "cont_test_ref";
		sim.init(1000, 1350);
		sim.simulate();
		sim.close();
	}

	// spinup run
	{
		Simulator sim(par_file);
		sim.continuePrevious = false;
		sim.expt_dir = "cont_test_spinup";
		sim.init(1000, 1200);
		sim.simulate();
		sim.close();
	}

	// continuation run
	{
		Simulator sim(par_file);
		sim.continuePrevious = true;
		sim.continueFrom_stateFile    = sim.parent_dir + "/" + "cont_test_spinup/pf_saved_state.txt"; 
		sim.continueFrom_configFile   = sim.parent_dir + "/" + "cont_test_spinup/pf_saved_config.ini"; 
		sim.expt_dir = "cont_test_main";
		sim.init(1000, 1350);
		sim.simulate();
		sim.close();
	}

	return 0;
}



// # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// # Input parameters for Simulations
// # If you change the order of parameters below, you will get what you deserve
// # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// > STRINGS 
// # # > DIR
// # Directories for data output
// # homeDir			/home/chethana/codes/Flare/tests		# home dir - no spaces allowed
// traitsFile      tests/data/Amz_trait_filled_HD.csv
// metFile         tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE.csv
// co2File         tests/data/CO2_ELE_AmzFACE2000_2100.csv

// outDir  		pspm_output2		    # output dir name
// exptName 	 	test_main_ref2							# expt name

// emgProps        AmzFACE_D_PFATE_ELE_HD.txt
// cwmAvg          AmzFACE_Y_mean_PFATE_ELE_HD.txt       
// cwmperSpecies   AmzFACE_Y_PFATE_ELE_HD.txt
// traits          traits_ELE_HD.txt

// solver          IEBT

// evolveTraits    yes

// saveState             yes
// savedStateFile        pf_saved_state.txt
// savedConfigFile       pf_saved_config.ini
// saveStateInterval     200

// continueFromState     null # pspm_output11/test_spinup/pf_saved_state.txt  # Set to null if fresh start desired
// continueFromConfig    null # pspm_output11/test_spinup/pf_saved_config.ini # Set to null if fresh start desired

// > SCALARS
// # ** 
// # ** Solver parameters
// # **
// resolution     5
// timestep       0.1
// T_cohort_insertion    1

// # **
// # ** Simulation parameters
// # **
// nSpecies       3
// year0          1000
// yearf          1350
// nPatches       3    

// # **
// # ** Core traits (default values)  
// # **
// lma            0.122 # 0.1163 # leaf mass per leaf area [kg/m2]
// zeta           0.20			# fine root mass per leaf area [kg/m2] 0.7
// fcr            0.47         # coarse root mass per unit stem mass
// hmat           29.18		# maximum height [m]
// fhmat          0.8          # height at reproductive maturity as fraction of hmat
// seed_mass      3.8e-5	    # [kg]
// wood_density   690	        # [kg/m3]
// p50_xylem      -2.29        # Xylem P50 [MPa]  -4.515

// # p50_leaf     -1.5		# Leaf hydraulic capacity [MPa]
// K_leaf         0.5e-16		# Leaf conductance [m]  (Tuned to gs)
// K_xylem        4e-16		# Leaf conductance [m]
// b_leaf         1			# Shape parameter of leaf vulnerabilty curve [-] (Ref: Joshi et al 2022)
// b_xylem        1            # Shape parameter of xylem vulnerabilty curve [-] 

// # **
// # ** Phydro paramaters  
// # **
// kphio          0.087       # Quantum yield efficiency
// alpha          0.095       # Cost of maintaining photosynthetic capacity (Ref: Joshi et al 2022, removed outliers Helianthus and Glycine)
// gamma          1.052       # Cost of maintaining hydraulic pathway  (Ref: Joshi et al 2022, removed outliers Helianthus and Glycine)     


// # **
// # ** Allocation and geometric paramaters  
// # **
// m   1.5	 # crown shape paramaters
// n   2
// fg  0.1  # upper canopy gap fraction

// a   75     # height-diameter allometry 114
// c   6000   # crown area allometry
// b   0	   # bark allometry


// # ** LAI model
// optimize_lai          0  # 1
// Cc                    0.3  # Leaf construction costs per unit mass 
// Chyd                  0.00
// response_intensity    3  # speed of LAI response. This is calibrated to give ~3 months response lag
// lai_deriv_step     1e-4  # stepsize to calculate profit derivative wrt LAI
// max_alloc_lai		0.5	 # max fraction of npp that can be allocated to LAI increment
// lai0                  1.8  # initial LAI

// # **
// # ** Leaf Economics
// # **
// # For all values below, Ref: Wang et al 2021 Leaf economics)
// les_u             768   # [dimensionless]
// les_k1             24.5   # g biomass / mol CO2 (see cbio below)
// les_k2           0.0864   # (mol-CO2/day) / (umol-CO2/s)
// les_hT_dH         65.33e3  # J mol-1
// les_hT_c          26.35   # - 
// les_molar_R       8.31    # J mol-1 K-1
// les_cc            13    # dimensionless    (Ref: Colin)

// # ** 
// # ** Respiration and turnover 
// # **
// rd  0.015                 # ratio of leaf dark respiration rate to vcmax [-]  (Ref: Colin)
// # rl  
// rr  0.6275                  # 5.3   # fine-root respiration rate [yr-1]  (kg/kg/yr) kg biomass lost = rr * kg in roots
// rs  0.032                # 0.16  # sapwood respiration rate [yr-1]    (kg/kg/yr) kg biomass lost = rs * kg in sapwood

// ll   0.5   # leaf lifespan [yr]
// lr   1     # fine root lifespan [yr]

// cbio 2.45e-2    # kg biomass per mol CO2
// y    0.75		# yield factor accounting for growth respiration Ref: Educated guess, also used in other models)

// k_light	0.5		# Light extinction coefficient

// # ** 
// # ** Demographics
// # **
// a_f1 0.15    # Max fraction of biomass invested into reproduction
// a_f2 10      # steepness of switch to reproduction

// ll_seed  15   # seed lifespan in dormant pool (seed avg dormancy)


// # **
// # ** Dispersal and germination
// # **
// Sd            1e-5     # probability of survival during dispersal
// npp_Sghalf    0.5      # 1.5 required productivity for 50% probability of survival during germination

// # **
// # ** Mortality
// # ** 
// # mI      0.002    # 0.002 natural mortality rate
// # mD      0.008    # 0.001 mortality rate due to diameter
// # mD_e    0.2     # exponent in mortality rate due to diameter
// # mS      1      # 0.0  scaling of rgr in mortality rate due to carbon starvation
// # mS0     1e-4   # mort = -log(mS0 + mS*rgr)

// # mS0    1    # used for mort = mS0*exp(-mS*rgr)
// # mS     10 # 5000      

// # Parameters of mortality function
// # exp(c0 + clnD*log(D) + cD*D + cWD*(wd*wd-cWD0*cWD0) + par.cS0*exp(-par.cS*bp.dmass_dt_tot))
// # c0     -5
// # clnD   -0.3
// # cD      0.2
// # cWD    -1.48
// # cWD0    0.690
// # cS0     0   # 1
// # cS      10

// # r = c0 + cL*log(L) + clnD*log(D) + cD*(D) + cG*log(rgr*D) + cWD*(WD - cWD0)
// # M = 1/(1+exp(-r))
//  c0     -5   # a
//  cL      0 # 0.15   # b
//  clnD    0 # -1   # c
//  cD      0.005  # d
//  cG      -2 # -0.35   # e
//  cWD    -0.05
//  cWD0    650
//  cS0     0

//  cD0     0.06  # cD0*D^2
//  cD1     0.05 # cD1*exp(-D/0.01)

//  m_alpha  0.0598
//  m_beta   18.7159
//  m_gamma  0.0094

// # **
// # ** Disturbance
// # **
// T_return              100       # return interval of disturbance
// T_seed_rain_avg       0.099       # years over which seed rain is averaged


// > ARRAYS
// # > PARAM_SWEEPS
// # parameter sets to loop over. These override any of the parameters set above.
// # all of these MUST BE SET AFTER particle-system initialization. 
// # 	  else, particleSystem will initialize with 1st value in these vectors by default 
// # c_full = 0.1 0.12 0.14 0.16 0.19 0.22 0.26 0.3 0.35 0.41 0.48 0.56 0.65 0.77 0.89 1.05 1.22 1.43 1.67 1.96 2.29 2.68 3.13 3.66 4.28 5 5.85 6.84 8 9.36 10.95 12.8  -1
// # cS_full = 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 -1
// # c = 0.1 0.14 0.19 0.26 0.35 0.48 0.65 0.89 1.22 1.67 2.29 3.13 4.28 5.85 8 10.95 -1
// # c_offset = 0.12 0.16 0.22 0.3 0.41 0.56 0.77 1.05 1.43 1.96 2.68 3.66 5 6.84 9.36 12.8 -1
// # bvec		0.0002 0.002 0.02 0.2 -1





