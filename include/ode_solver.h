#ifndef PSPM_ODE_SOLVER_H_
#define PSPM_ODE_SOLVER_H_

#include "rkck45.h"
#include "lsoda.h"

#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <fstream>

enum SolverType {ODE_RKCK45, ODE_LSODA};

class OdeSolver{
	private:
	SolverType type;
	void * solver = nullptr;
	int nfe_cumm = 0;

	public:
	struct{
		double abs_tol = 1e-8;
		double rel_tol = 1e-8;
	} control;

	private:
	inline void deleteSolver(){
		if      (type == ODE_RKCK45) delete static_cast<RKCK45*>(solver);	
		else if (type == ODE_LSODA)  delete static_cast<LSODA*>(solver);
		solver = nullptr;	
	}

	inline void createSolver(double t_start, double rtol, double atol){
		if      (type == ODE_RKCK45) solver = new RKCK45(t_start, rtol, 1e-8);
		else if (type == ODE_LSODA)  solver = new LSODA();
	}

	public:
	OdeSolver(std::string method, double t_start, double rtol, double atol){
		if      (method == "rk45ck") type = ODE_RKCK45;
		else if (method == "lsoda")  type = ODE_LSODA;
		else throw std::runtime_error("Fatal: Unknown ODE method " + method);
		
		createSolver(t_start, rtol, atol);
	}

	~OdeSolver(){
		deleteSolver();
	}

	// FIXME: Implement copy-constructor

	// TODO: Implement this via copy & swap idiom
	OdeSolver& operator=(const OdeSolver &rhs){
		if (this != &rhs){
			std::cout << "OdeSolver::operator= entered\n"; 
			// copy construct the new solver as per the rhs's type
			// do this first, so that if this throws an exception, the lhs object is still valid
			void * solver2;
			std::cout << "RKCK45 constructor entered: " << solver << '\n';
			if      (rhs.type == ODE_RKCK45) solver2 = new RKCK45(*static_cast<RKCK45*>(rhs.solver));
			else if (rhs.type == ODE_LSODA)  solver2 = new LSODA(*static_cast<LSODA*>(rhs.solver));
			else throw std::runtime_error("Fatal: Unknown ODE Solver type " + type);

			this->deleteSolver(); // delete current solver (based on current type)

			// copy new solver type and other meta from rhs
			type = rhs.type;
			control = rhs.control;
			nfe_cumm = rhs.nfe_cumm;
			solver = solver2;
		}
		return *this;
	}

	void reset(double t_start, double rtol, double atol){
		nfe_cumm = 0;
		control.abs_tol = atol;
		control.rel_tol = rtol;
		deleteSolver();
		createSolver(t_start, rtol, atol);
	}


	template <class Functor, class AfterStep>
	void step_to(double t_stop, double &t, std::vector<double>&y, Functor &derivs, AfterStep &after_step){
		if (t_stop == t || y.size() == 0) return;
		
		if (type == ODE_RKCK45){
			RKCK45 * sol = static_cast<RKCK45*>(solver);
			sol->Step_to(t_stop, t, y, derivs, after_step);
			nfe_cumm += sol->get_fn_evals();
		}
		else if (type == ODE_LSODA ){
			LSODA* sol = static_cast<LSODA*>(solver);
			sol->set_istate(1); // forces re-initialization
			sol->lsoda_update(derivs, after_step, y.size(), y, &t, t_stop, nullptr, control.rel_tol, control.abs_tol);
			if(sol->get_istate() <= 0) {
				std::cerr << "LSODA Error: istate = " << sol->get_istate() << std::endl;
				return;
			}
			nfe_cumm += sol->get_fncalls();
		}
		
		//else if (type == ODE_LSODA ){
			////LSODA* sol = (LSODA*)solver;
			////delete sol;
			//LSODA * sol = new LSODA();
			//sol->lsoda_update(derivs, y.size(), y, &t, t_stop, nullptr, control.rel_tol, control.abs_tol);
			//if(sol->get_istate() <= 0) {
				//std::cerr << "LSODA Error: istate = " << sol->get_istate() << std::endl;
				//return;
			//}
			//nfe_cumm += sol->get_fncalls();
			//delete sol;
		//}
	}

	int get_fn_evals(){
		if      (type == ODE_RKCK45) return static_cast<RKCK45*>(solver)->get_fn_evals();	
		else if (type == ODE_LSODA)  return nfe_cumm;
		else return -1;
	}

	void save(std::ostream &fout){
		fout << "odeSolver::v1\n";

		fout << control.abs_tol << ' '
		     << control.rel_tol << ' ';

		fout << static_cast<int>(type) << '\n';

		if      (type == ODE_RKCK45) static_cast<RKCK45*>(solver)->save(fout);	
		else if (type == ODE_LSODA)  throw std::runtime_error("Cannot save the state for LSODA solver.");
	}

	void restore(std::istream &fin){
		deleteSolver(); // delete current solver as its type may be different from the saved type...

		std::string s; fin >> s; // discard version number

		fin >> control.abs_tol
		    >> control.rel_tol;

		int m;
		fin >> m; 
		type = SolverType(m);

		// ...then recreate solver based on saved type
		createSolver(0,0,0); // dummy arguments here are fine, as all variables will be recreated from saved file

		if      (type == ODE_RKCK45) static_cast<RKCK45*>(solver)->restore(fin);	
		else if (type == ODE_LSODA)  throw std::runtime_error("Cannot restore the state for LSODA solver at this point.");
	}

	void printODEsolvermethod(){
		std::cout << "ODE solver: " << type << std::endl;
	}


};



#endif

