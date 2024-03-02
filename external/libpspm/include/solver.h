#ifndef  PSPM_PSPM_SOLVER_H_
#define  PSPM_PSPM_SOLVER_H_

#include <vector>
#include <list>
#include <map>
#include <string>
#include <random>

#include "environment_base.h"
#include "species.h"
#include "ode_solver.h"
#include "cubic_spline.h"

enum PSPM_SolverType {SOLVER_FMU, 
                      SOLVER_MMU, 
                      SOLVER_CM, 
                      SOLVER_EBT, 
                      SOLVER_IFMU, 
                      SOLVER_ABM, 
                      SOLVER_IEBT,
                      SOLVER_ICM};

class Solver{
	// define the type of a pointer to member functions calcRates_XXX()
	using calcRatesPointer = void (Solver::*)(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	// define the type of a pointer to member functions stepU_XXX()
	using stepUPointer = void (Solver::*)(double t, std::vector<double>& S, std::vector<double>& dSdt, double dt);

	private:
	static std::map<std::string, PSPM_SolverType> methods_map;

	bool debug = false;
	PSPM_SolverType method;

	// int n_statevars_internal = 0;		// number of internal i-state variables (x and/or u)

	std::default_random_engine generator; // random number generator (used by ABM)

	public:	
	OdeSolver odeStepper;
	EnvironmentBase * env = nullptr;
	int n_statevars_system = 0;			// number of s-state variables 
	std::vector<double> s_state;
	
	// The current state of the system, {t, S, dS/dt} 
	// FIXME: Setting initial value of current_time to anything other than 0 breaks Falster17 seed rain in Plant model. Investigate. 
	//  ^ Solution: init of current_time must happen before species_initialization. Maybe via an argument to solver constructor?
	double current_time = 0;        // these are synced with ODE solver only after a successful step
	std::vector <double> state;     // +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 

	std::vector<Species_Base*> species_vec;	

	// FIXME: Should these control params be saved when saving solver state?
	struct{
		double ode_eps = 1e-6;
		double ode_initial_step_size = 1e-6;
		double convergence_eps = 1e-6;
		std::vector<double> cm_grad_dx = {1e-6};
		bool update_cohorts = true;
		bool cm_remove_cohorts = true;
		int  max_cohorts = 500;
		double cm_dxcut = 1e-10;
		double ebt_ucut = 1e-10;
		double ebt_grad_dx = 1e-6;
		double ebt_merge_dxcut = 0; //1e-6;
		double ode_rk4_stepsize = 0.1;
		double ode_ifmu_stepsize = 0.1;
		bool ifmu_centered_grids = true;
		bool integral_interpolate = true;
		double ifmu_order = 1;
		double abm_n0 = 100;
		double abm_burnin = 1000;
		int abm_numChains = 4;
		double abm_stepsize = 0.02;
		bool abm_init_on_grid = true;
		double cohort_insertion_dt = 1e20;
		double cohort_insertion_tol = 1e-12;
		bool sync_cohort_insertion = false;
		bool cm_use_log_densities = true;
	} control;
	

	private:
	std::vector<double> pi0;
	double N0;
	double t_next_cohort_insertion = 0;

	void realizeEbtBoundaryCohort(Species_Base * spp);
	// void realizeEbtnBoundaryCohort(Species_Base * spp);
	void restoreEbtBoundaryCohort(Species_Base * spp);
	// void restoreEbtnBoundaryCohort(Species_Base * spp);

	public:	
	Solver(PSPM_SolverType _method, std::string ode_method = "rk45ck");
	Solver(std::string _method, std::string ode_method = "rk45ck");

	void addSystemVariables(const std::vector<double>& s_state_0);
	void addSpecies(std::vector<int> _J, std::vector<double> _xb, std::vector<double> _xm, std::vector<bool> log_breaks, Species_Base* _mod, int _n_accumulators, double input_birth_flux = -1);
	void addSpecies(std::vector<std::vector<double>> xbreaks, Species_Base* _mod, int _n_accumulators, double input_birth_flux = -1);
	void removeSpecies(Species_Base* spp);

	//Species<Model>* get_species(int id);
	int n_species();
	int n_statevars(Species_Base * spp);
	int n_statevars_cohort(Species_Base * spp);

	void setEnvironment(EnvironmentBase * _env);

	// void resetState(double t0 = 0);
	void initialize(double t0 = 0);
	void resizeStateFromSpecies();

	void initializeSpecies(Species_Base * s);

	void copyStateToCohorts(std::vector<double>::iterator state_begin);		////const int size();
	void copyCohortsToState();
	
	std::vector<double> maxState(int species_id);

	void updateEnv(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);

	/// @brief calculate \f$du/dt\f$ using the FMU solver
	void calcRates_FMU(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	
	void calcOdeRatesImplicit(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	void stepAccumulators(double dt); // this will step accumulators only to current_tile+dt. This is used by implicit solvers
	void stepSystemVars(const std::vector<double> &sys_rates_prev, double dt);

	void stepU_iFMU(double t, std::vector<double> &S, std::vector<double> &dSdt, double dt);
	void stepU_iEBT(double t, std::vector<double> &S, std::vector<double> &dSdt, double dt);
	void stepU_iCM(double t, std::vector<double> &S, std::vector<double> &dSdt, double dt);

	/// @brief calculate \f$du/dt\f$ using the EBT solver
	void calcRates_EBT(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	void addCohort_EBT();
	void removeDeadCohorts_EBT();
	void mergeCohorts_EBT();

	/// @brief calculate \f$du/dt\f$ using the CM solver
	void calcRates_CM(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	void addCohort_CM();
	void removeCohort_CM();

	void stepABM(double t, double dt);
	
	template<typename AfterStepFunc>
	void stepTo_explicit(double tstop, AfterStepFunc & afterStep, calcRatesPointer rates_func);

	template<typename AfterStepFunc>
	void stepTo_implicit(double tstop, AfterStepFunc &afterStep, stepUPointer step_u_func);

	template<typename AfterStepFunc>
	void stepTo_abm(double tstop, AfterStepFunc & afterStep);

	template<typename AfterStepFunc>
	void step_to(double tstop, AfterStepFunc &afterStep_user);

	void step_to(double tstop);

	/// @brief Calculate total flux of newborns, without updating environment
	/// This function is used internally by the solver, because the solver explicitly manages environment computation. 
	/// In the solver, all E,g,m,f,B are all calculated at the beginning of the step, in that order, 
	/// so we dont want to again update E during the calculation of B
	double calcSpeciesBirthFlux(int k, double t); // TODO: Make this private. Because this does not update Env, so unsafe for users to call. Users should call newborns_out() instead because it does updateEnv + calcSpeciesBirthFlux()
	
	/// @brief Calculate the total flux of newborns given the current population structure, optionally after recomputing environment  
	/// @param t             current time
	/// @param recompute_env whether environment should be recomputed before calculating fecundity rates
	/// @return              birth flux
	/// This function is the same as calcSpeciesBirthFlux but provides the option to update env before computing f,B. 
	/// Typically, this function is meant to be used at the end of a timestep, like so: E,g,m,f,B,u -----> u'--> afterStep: {E',B'}
	/// For consistency, i.e., B' ~ u',E' (not ~ u',E), it's necessary to update env. 
	std::vector<double> newborns_out(double t, bool recompute_env = true);  // This is the actual system reproduction (fitness) hence biologically relevant
	std::vector<double> u0_out(double t);        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value

	void print();
	
	// // integrals over size and density
	// /// @brief Computes the integral \f[I = \int_{x_b}^{x_m} w(z,t)u(z)dz\f] for the specified species. For details, see @ref integrate_wudx_above
	// template<typename wFunc>
	// double integrate_x(wFunc w, double t, int species_id);

	/// @brief Computes the partial integral \f[I = \int_{x_{low}}^{x_m} w(z,t)u(z)dz\f] for the specified species. 
	template<typename wFunc>
	double integrate_wudx_above(wFunc w, double t, const std::vector<double>& xlow, int species_id);

	// template<typename wFunc>
	// double integrate_wudxn_above(wFunc w, double t, std::vector<double> xnlow, std::vector<double> xnhigh, int species_id);

	template<typename wFunc>
	double state_integral(wFunc w, double t, int species_id);

	std::vector<double> getDensitySpecies1D(int k, int dim, const std::vector<double>& breaks, Spline::Extr extrapolation_method = Spline::ZERO);
	// std::vector<std::vector<double>> getDensitySpecies2D(int k, const std::vector<int>& axes, const std::vector<std::vector<double>>& breaks, Spline::Extr extrapolation_method);
	
	void save(std::ostream &fout);
	void restore(std::istream &fin, std::vector<Species_Base*> spp_proto);

};

#include "../src/solver.tpp"
#include "../src/size_integrals.tpp"


#endif



