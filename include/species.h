#ifndef  PSPM_PSPM_SPECIES_H_
#define  PSPM_PSPM_SPECIES_H_

#include <vector>
#include <list>
#include <string>

//#include "iterator_set.h"
#include "cohort.h"

// forward declaration of Solver so Species can befriend it
class Solver;

class Species_Base{
	// Solver should be able to access Species' privates
	friend class Solver;

	protected: // private members
	int J = 0;	

	// std::list<double> birth_flux_out_history;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	
	// These are only used by FMU solver. 
	// Others derive them from the state. 
	// Kept private so users dont accidently access them for other solvers
	std::vector<int> dim_centres, dim_edges;
	int n_grid_centres, n_grid_edges; // These are cummulative products of dim_centres and dim_edges respectively
	std::vector <std::vector<double>> X; // multiple state model
	std::vector <std::vector<double>> x; // multiple state model
	std::vector <std::vector<double>> h;
	// std::vector <double> schedule; // used only by CM/EBT

	double noff_abm = 0; // used by ABM solver to insert offspring

	public:
	size_t istate_size;
	int n_accumulators = 0;
	double birth_flux_in = -9.9e20; //std::numeric_limits<double>::quiet_NaN();
	
	//debug only
	bool bfin_is_u0in = false;
	
	std::vector <double> xb;
	
	public: // public functions
	virtual ~Species_Base() = 0;
	
	// int cohortsize(); // 
	
	int xsize();
	int stateSizeTotal();

	void set_inputBirthFlux(double b);
	void set_bfin_is_u0in(bool flag);


	public:
	virtual void clear_vectors() = 0;
	virtual void resize(int _J) = 0;
	virtual std::vector<double> get_maxSize(int skip) = 0;
	virtual void print() = 0;
	virtual void print_extra(); // not pure virtual, by defualt, there is nothing extra to print.

	virtual void set_xb(std::vector<double> _xb) = 0;
	virtual void set_ub(double _ub) = 0;
	virtual void set_birthTime(int i, double t0) = 0;
	virtual void setX(int i, std::vector<double> _x) = 0;
	virtual void setU(int i, double _u) = 0;

	virtual std::vector<double> getX(int i) = 0;
	virtual double getU(int i) = 0;

	// virtual double dXn (int i) = 0;
	// virtual double dXn(std::vector<double> xn1, std::vector<double> xn2) = 0;
	// virtual std::vector<double> cohort_dist(std::vector<double> xn1, std::vector<double> xn2) = 0;
	// virtual double next_xk_desc(double xkn, int k) = 0;
	// virtual double next_xk_asc(double xkn, int k) = 0;
	// virtual std::vector<double> next_xn_desc(std::vector<double> xn) = 0;
	// virtual std::vector<double> next_xn_asc(std::vector<double> xn) = 0;

	virtual void initBoundaryCohort(double t, void * env) = 0;
	virtual double init_density(int i, void * env) = 0; // FIXME: Should init_density take t as input?

	virtual void initAccumulators(double t, void * env) = 0;
	virtual void initAndCopyAccumulators(double t, void * env, std::vector<double>::iterator &it) = 0;
	virtual void copyAccumulatorsToCohorts(std::vector<double>::iterator &it) = 0;
	virtual void copyAccumulatorsToState(std::vector<double>::iterator &it) = 0;
	virtual void accumulatorRates(std::vector<double>::iterator &it) = 0;

	virtual double establishmentProbability(double t, void * env) = 0;
	virtual double calc_boundary_u(std::vector<double> gb, double pe) = 0;
	virtual double get_boundary_u() = 0;

	virtual void triggerPreCompute() = 0;

	// // TODO: argument x can probably be removed from these functions
	virtual std::vector<double> growthRate(int i, double t, void * env) = 0;
	virtual std::vector<double> growthRateOffset(int i, const std::vector<double>& x, double t, void * env) = 0;
	virtual std::vector<std::vector<double>> growthRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx) = 0;
	// // virtual std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env) = 0;
	virtual double mortalityRate(int i, double t, void * env) = 0;
	virtual std::vector<double> mortalityRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx) = 0;
	virtual double birthRate(int i, double t, void * env) = 0;

	virtual void addCohort(int n = 1) = 0;
	template<class T> void addCohort(T bc);

	virtual void markCohortForRemoval(int i) = 0;
	virtual void markDensestCohort() = 0;
	virtual void markDenseCohorts(double dxcut) = 0;
	virtual void markDeadCohorts(double ucut) = 0;
	virtual void removeMarkedCohorts() = 0;
	// virtual void mergeCohortsAddU(std::vector<double> dxcut) = 0;

	virtual void sortCohortsDescending(size_t dim, int skip=0) = 0;
	virtual void sortCohortsAscending(size_t dim, int skip=0) = 0;
	
	virtual void save(std::ostream &fout) = 0;
	virtual void restore(std::istream &fin) = 0;

	virtual void printCohortVector(std::ostream &out) = 0;

//	virtual void backupCohort(int j) = 0;
//	virtual void restoreCohort(int j) = 0;
//	virtual void copyBoundaryCohortTo(int j) = 0;
};



template <class Model>
class Species : public Species_Base{
	protected:
	std::vector<Cohort<Model>> cohorts;
	Cohort<Model> boundaryCohort;
	
	//Cohort<Model> savedCohort; // a cohort to save a backup of any other cohort

	public:
	// Make these virtual? - not needed. They are virtual by default.
	// Species(std::vector<double> breaks = std::vector<double>());
	
	// Constructor based on model instance
	// Species(const Model& M);
	
	// Constructor based on arguments required by model constructor
	template <typename... ARGS>
	Species(ARGS... args);

	// Species<Model>* create(); // virtual constructor needed for deserialization

	
	void clear_vectors();
	void resize(int _J);
	std::vector<double> get_maxSize(int skip=0);
	void print();
	using Species_Base::print_extra;

	void set_xb(std::vector<double> _xb);
	void set_ub(double _ub);
	void set_birthTime(int i, double t0);
	void setX(int i, std::vector<double> _x);
	void setU(int i, double _u);

	std::vector<double> getX(int i);
	double getU(int i);

	// double dXn (int i);
	// double dXn(std::vector<double> xn1, std::vector<double> xn2);
	// std::vector<double> cohort_dist(std::vector<double> xn1, std::vector<double> xn2);
	// double next_xk_desc(double xkn, int k);
	// double next_xk_asc(double xkn, int k);
	// std::vector<double> next_xn_desc(std::vector<double> xn);
	// std::vector<double> next_xn_asc(std::vector<double> xn);

	void initBoundaryCohort(double t, void * env);
	double init_density(int i, void * env);

	void initAccumulators(double t, void * env);
	void initAndCopyAccumulators(double t, void * env, std::vector<double>::iterator &it);
	void copyAccumulatorsToCohorts(std::vector<double>::iterator &it);
	void copyAccumulatorsToState(std::vector<double>::iterator &it);
	void accumulatorRates(std::vector<double>::iterator &it);
	
	double establishmentProbability(double t, void * env);
	double calc_boundary_u(std::vector<double> gb, double pe);
	double get_boundary_u();

	void triggerPreCompute();

	// // TODO: argument x can probably be removed from these functions
	std::vector<double> growthRate(int i, double t, void * env);
	std::vector<double> growthRateOffset(int i, const std::vector<double>& x, double t, void * env);
	std::vector<std::vector<double>> growthRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx);
	// // std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env);
	double mortalityRate(int i, double t, void * env);
	std::vector<double> mortalityRateGradient(int i, double t, void * env, const std::vector<double>& grad_dx);
	double birthRate(int i, double t, void * env);

	void addCohort(int n = 1);
	template<class T> void addCohort(T bc);

	void markCohortForRemoval(int i);
	// JJ Note: A sorting-independent way to implement the following would be to use FoF algo. 
	//          Postponing this to later. For now, these cohorts do nothing. In any case, they are only required by CM, 
	//          and for current purposes, its ok if cohorts are not removed
	void markDensestCohort();
	void markDenseCohorts(double dxcut);
	void markDeadCohorts(double ucut);
	void removeMarkedCohorts();
	// void mergeCohortsAddU(std::vector<double> dxcut);

	void sortCohortsDescending(size_t dim, int skip=0);
	void sortCohortsAscending(size_t dim, int skip=0);
	
	void save(std::ostream &fout);
	void restore(std::istream &fin);

	void printCohortVector(std::ostream &out);

//	void backupCohort(int j);
//	void restoreCohort(int j);
//	void copyBoundaryCohortTo(int j);

	public:
	Cohort<Model>& getCohort(int i);
};


#include "../src/species.tpp"

#endif
