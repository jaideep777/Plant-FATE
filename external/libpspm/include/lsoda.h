/***
 *        Created:  2018-08-14

 *         Author:  Dilawar Singh <dilawars@ncbs.res.in>
 *   Organization:  NCBS Bangalore
 *        License:  MIT License
 *
 *       Modified:  2021-10-05
 *         Author:  Jaideep Joshi <jaideep777@gmail.com>
 *   Organization:  IIASA Vienna
 */

#ifndef PSPM_ODE_LSODA_H_
#define PSPM_ODE_LSODA_H_

#include <array>
#include <cmath>
#include <memory>
#include <vector>

//using namespace std;

/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Type definition of LSODA ode system. See the file test_LSODA.cpp
 * for an example.
 *
 * @Param time, double
 * @Param y, array of double.
 * @Param dydt, array of double
 * @Param data, void*
 *
 * @Returns void
 */
/* ----------------------------------------------------------------------------*/
//typedef void (*LSODA_ODE_SYSTEM_TYPE)(double t, double *y, double *dydt, void *);

class LSODA {

public:
    LSODA();
    ~LSODA();

    size_t idamax1(const std::vector<double> &dx, const size_t n, const size_t offset);

    void dscal1(const double da, std::vector<double> &dx, const size_t n, const size_t offset);

    double ddot1(const std::vector<double> &a, const std::vector<double> &b, const size_t n,
        const size_t offsetA, const size_t offsetB);

    void daxpy1(const double da, const std::vector<double> &dx, std::vector<double> &dy,
        const size_t n, const size_t offsetX, const size_t offsetY);

    void dgesl(const std::vector<std::vector<double>> &a, const size_t n, std::vector<int> &ipvt,
        std::vector<double> &b, const size_t job);

    void dgefa(
        std::vector<std::vector<double>> &a, const size_t n, std::vector<int> &ipvt, size_t *const info);

    template <class Functor>
	void prja(const size_t neq, std::vector<double> &y, Functor &derivs, void *_data);

    template <class Functor, class AfterStep>
    void lsoda(Functor &derivs, AfterStep &after_step, const size_t neq, std::vector<double> &y, double *t,
        double tout, int itask, int *istate, int iopt, int jt, std::array<int, 7> &iworks,
        std::array<double, 4> &rworks, void *_data);

    template <class Functor>
    void correction(const size_t neq, std::vector<double> &y, Functor &derivs,
        size_t *corflag, double pnorm, double *del, double *delp, double *told,
        size_t *ncf, double *rh, size_t *m, void *_data);

    template <class Functor>
    void stoda(const size_t neq, std::vector<double> &y, Functor &derivs, void *_data);

    // We call this function in VoxelPools::
    template <class Functor, class AfterStep>
    void lsoda_update(Functor &derivs, AfterStep &after_step, const size_t neq, std::vector<double> &y,
        double *t, const double tout, 
        void *const _data, double rtol = 1e-6, double atol = 1e-6  // Tolerance
    );

    template<class AfterStep>
    void successreturn(AfterStep &after_step,
        std::vector<double> &y, double *t, int itask, int ihit, double tcrit, int *istate);

    void terminate(int *istate);
    void terminate2(std::vector<double> &y, double *t);
    void _freevectors(void);
    void ewset(const std::vector<double> &ycur);
    void resetcoeff(void);
    void solsy(std::vector<double> &y);
    void endstoda(void);
    void orderswitch(
        double *rhup, double dsm, double *pdh, double *rh, size_t *orderflag);
    void intdy(double t, int k, std::vector<double> &dky, int *iflag);
    void corfailure(double *told, double *rh, size_t *ncf, size_t *corflag);
    void methodswitch(double dsm, double pnorm, double *pdh, double *rh);
    void cfode(int meth_);
    void scaleh(double *rh, double *pdh);
    double fnorm(int n, const std::vector<std::vector<double>> &a, const std::vector<double> &w);
    double vmnorm(const size_t n, const std::vector<double> &v, const std::vector<double> &w);

    static bool abs_compare(double a, double b);

	void resize_system(int nyh, int lenyh);	
	void set_istate(int n);
	int get_istate() const;
	int get_fncalls();

private:
	double ETA = 2.2204460492503131e-16;
    
	size_t ml, mu, imxer;
    double sqrteta;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here. With newer compiler,
    // initialization can be done here.
    std::array<size_t, 3> mord;
    std::array<double, 13> sm1;

    std::array<double, 14> el;   // = {0};
    std::array<double, 13> cm1;  // = {0};
    std::array<double, 6> cm2;   // = {0};

    std::array<std::array<double, 14>, 13> elco;
    std::array<std::array<double, 4>, 13> tesco;

    size_t illin, init, ierpj, iersl, jcur, l, miter, maxord, maxcor, msbp, mxncf;

    int kflag, jstart;

    size_t ixpr = 0, jtyp, mused, mxordn, mxords = 12;
    size_t meth_;

    size_t n, nq, nst, nfe, nje, nqu;
    size_t mxstep, mxhnil;
    size_t nslast, nhnil, ntrep, nyh;

    double ccmax, el0, h_ = .0;
    double hmin, hmxi, hu, rc, tn_ = 0.0;
    double tsw, pdnorm;
    double conit, crate, hold, rmax;

    size_t ialth, ipup, lmax;
    size_t nslp;
    double pdest, pdlast, ratio;
    int icount, irflag;

    std::vector<double> ewt;
    std::vector<double> savf;
    std::vector<double> acor;
    std::vector<std::vector<double>> yh_;
    std::vector<std::vector<double>> wm_;

    std::vector<int> ipvt;

private:
    int itol_ = 2;
    std::vector<double> rtol_;
    std::vector<double> atol_;

private: // convenience members so that LSODA class is entirely standalone
	int iState = 1;  // LSODA should always be started with istate=1
	//int nEq;		 // number of ODEs
	std::vector<double> yout;

public:
    void *param = nullptr;
};


#include "../src/lsoda.tpp"

#endif /* end of include guard: LSODE_H */
