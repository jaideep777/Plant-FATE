/*
 * HISTORY:
 * This is a CPP version of the LSODA library for integration into MOOSE
 somulator.
 * The original was aquired from
 * http://www.ccl.net/cca/software/SOURCES/C/kinetics2/index.shtml and modified
 by
 * Heng Li <lh3lh3@gmail.com>. Heng merged several C files into one and added a
 * simpler interface. [Available
 here](http://lh3lh3.users.sourceforge.net/download/lsoda.c)

 * The original source code came with no license or copyright
 * information. Heng Li released his modification under the MIT/X11 license. I
 * maintain the same license. I have removed quite a lot of text/comments from
 * this library. Please refer to the standard documentation.
 *
 * Contact: Dilawar Singh <dilawars@ncbs.res.in>
*/

#include "lsoda.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

//#include "helper.h"

using namespace std;


LSODA::LSODA()
{
	// Initialize arrays.
	mord = {{12, 5}};
	sm1  = {{0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025}};
	el   = {{0}};
	cm1  = {{0}};
	cm2  = {{0}};
}

LSODA::~LSODA()
{
}

bool LSODA::abs_compare(double a, double b)
{
	return (std::abs(a) < std::abs(b));
}

/* Purpose : Find largest component of double vector dx */
size_t LSODA::idamax1(const vector<double> &dx, const size_t n, const size_t offset = 0)
{

	size_t v = 0, vmax = 0;
	size_t idmax = 1;
	for(size_t i = 1; i <= n; i++) {
		v = abs(dx[i + offset]);
		if(v > vmax) {
			vmax  = v;
			idmax = i;
		}
	}
	return idmax;

	// Following has failed with seg-fault. Probably issue with STL.
	// return std::max_element( dx.begin()+1+offset, dx.begin()+1+n, LSODA::abs_compare) -
	// dx.begin() - offset;
}

/* Purpose : scalar vector multiplication
   dx = da * dx
*/
void LSODA::dscal1(
	const double da, vector<double> &dx, const size_t n, const size_t offset = 0)
{
	// FIXME: n is not used here. why?
	(void)n;

	std::transform(dx.begin() + 1 + offset, dx.end(), dx.begin() + 1 + offset,
		[&da](double x) -> double { return da * x; });
}

/* Purpose : Inner product dx . dy */
double LSODA::ddot1(const vector<double> &a, const vector<double> &b, const size_t n,
	const size_t offsetA = 0, const size_t offsetB = 0)
{
	double sum = 0.0;
	for(size_t i = 1; i <= n; i++)
		sum += a[i + offsetA] * b[i + offsetB];
	return sum;
}

void LSODA::daxpy1(const double da, const vector<double> &dx, vector<double> &dy,
	const size_t n, const size_t offsetX = 0, const size_t offsetY = 0)
{

	for(size_t i = 1; i <= n; i++)
		dy[i + offsetY] = da * dx[i + offsetX] + dy[i + offsetY];
}

// See BLAS documentation. The first argument has been changed to vector.
void LSODA::dgesl(const vector<vector<double>> &a, const size_t n, vector<int> &ipvt,
	vector<double> &b, const size_t job)
{
	size_t k, j;
	double t;

	/*
	   Job = 0, solve a * x = b.
	*/
	if(job == 0) {
		/*
		   First solve L * y = b.
		*/
		for(k = 1; k <= n; k++) {
			t    = ddot1(a[k], b, k - 1);
			b[k] = (b[k] - t) / a[k][k];
		}
		/*
		   Now solve U * x = y.
		*/
		for(k = n - 1; k >= 1; k--) {
			b[k] = b[k] + ddot1(a[k], b, n - k, k, k);
			j    = ipvt[k];
			if(j != k) {
				t    = b[j];
				b[j] = b[k];
				b[k] = t;
			}
		}
		return;
	}
	/*
	   Job = nonzero, solve Transpose(a) * x = b.

	   First solve Transpose(U) * y = b.
	*/
	for(k = 1; k <= n - 1; k++) {
		j = ipvt[k];
		t = b[j];
		if(j != k) {
			b[j] = b[k];
			b[k] = t;
		}
		daxpy1(t, a[k], b, n - k, k, k);
	}
	/*
	   Now solve Transpose(L) * x = y.
	*/
	for(k = n; k >= 1; k--) {
		b[k] = b[k] / a[k][k];
		t    = -b[k];
		daxpy1(t, a[k], b, k - 1);
	}
}

// See BLAS documentation. All double* has been changed to std::vector .
void LSODA::dgefa(
	vector<vector<double>> &a, const size_t n, vector<int> &ipvt, size_t *const info)
{
	size_t j = 0, k = 0, i = 0;
	double t = 0.0;

	/* Gaussian elimination with partial pivoting.   */

	*info = 0;
	for(k = 1; k <= n - 1; k++) {
		/*
		   Find j = pivot index.  Note that a[k]+k-1 is the address of
		   the 0-th element of the row vector whose 1st element is a[k][k].
		*/
		j       = idamax1(a[k], n - k + 1, k - 1) + k - 1;
		ipvt[k] = j;
		/*
		   Zero pivot implies this row already triangularized.
		*/
		if(a[k][j] == 0.) {
			*info = k;
			continue;
		}
		/*
		   Interchange if necessary.
		*/
		if(j != k) {
			t       = a[k][j];
			a[k][j] = a[k][k];
			a[k][k] = t;
		}
		/*
		   Compute multipliers.
		*/
		t = -1. / a[k][k];
		dscal1(t, a[k], n - k, k);

		/*
		   Column elimination with row indexing.
		*/
		for(i = k + 1; i <= n; i++) {
			t = a[i][j];
			if(j != k) {
				a[i][j] = a[i][k];
				a[i][k] = t;
			}
			daxpy1(t, a[k], a[i], n - k, k, k);
		}
	} /* end k-loop  */

	ipvt[n] = n;
	if(a[n][n] == 0.)
		*info = n;
}

/* Terminate lsoda due to illegal input. */
void LSODA::terminate(int *istate)
{
	if(illin == 5)
		cerr << "[lsoda] repeated occurrence of illegal input. run aborted.. "
				"apparent infinite loop."
			 << endl;
	else {
		illin++;
		*istate = -3;
	}
}

/* Terminate lsoda due to various error conditions. */
void LSODA::terminate2(vector<double> &y, double *t)
{
	for(size_t i = 1; i <= n; i++)
		y[i] = yh_[1][i];
	*t    = tn_;
	illin = 0;
	return;
}



void LSODA::ewset(const vector<double> &ycur)
{
	switch(itol_) {
		case 1:
			for(size_t i = 1; i <= n; i++)
				ewt[i] = rtol_[1] * fabs(ycur[i]) + atol_[1];
			break;
		case 2:
			for(size_t i = 1; i <= n; i++)
				ewt[i] = rtol_[1] * fabs(ycur[i]) + atol_[i];
			break;
		case 3:
			for(size_t i = 1; i <= n; i++)
				ewt[i] = rtol_[i] * fabs(ycur[i]) + atol_[1];
			break;
		case 4:
			for(size_t i = 1; i <= n; i++)
				ewt[i] = rtol_[i] * fabs(ycur[i]) + atol_[i];
			break;
	}

} /* end ewset   */

/*
   Intdy computes interpolated values of the k-th derivative of the
   dependent variable vector y, and stores it in dky.  This routine
   is called within the package with k = 0 and *t = tout, but may
   also be called by the user for any k up to the current order.
   ( See detailed instructions in the usage documentation. )

   The computed values in dky are gotten by interpolation using the
   Nordsieck history array yh_.  This array corresponds uniquely to a
   vector-valued polynomial of degree nqcur or less, and dky is set
   to the k-th derivative of this polynomial at t.
   The formula for dky is

			 q
   dky[i] = sum c[k][j] * ( t - tn_ )^(j-k) * h_^(-j) * yh_[j+1][i]
			j=k

   where c[k][j] = j*(j-1)*...*(j-k+1), q = nqcur, tn_ = tcur, h_ = hcur.
   The quantities nq = nqcur, l = nq+1, n = neq, tn_, and h_ are declared
   static globally.  The above sum is done in reverse order.
   *iflag is returned negative if either k or t is out of bounds.
*/
void LSODA::intdy(double t, int k, vector<double> &dky, int *iflag)
{
	int ic, jp1 = 0;
	double c, r, s, tp;

	*iflag = 0;
	if(k < 0 || k > (int)nq) {
		fprintf(stderr, "[intdy] k = %d illegal\n", k);
		*iflag = -1;
		return;
	}
	tp = tn_ - hu - 100. * ETA * (tn_ + hu);
	if((t - tp) * (t - tn_) > 0.) {
		fprintf(
			stderr, "intdy -- t = %g illegal. t not in interval tcur - hu to tcur\n", t);
		*iflag = -2;
		return;
	}
	s  = (t - tn_) / h_;
	ic = 1;
	for(size_t jj = l - k; jj <= nq; jj++)
		ic *= jj;
	c = (double)ic;
	for(size_t i = 1; i <= n; i++)
		dky[i] = c * yh_[l][i];

	for(int j = nq - 1; j >= k; j--) {
		jp1 = j + 1;
		ic  = 1;
		for(int jj = jp1 - k; jj <= j; jj++)
			ic *= jj;
		c = (double)ic;

		for(size_t i = 1; i <= n; i++)
			dky[i] = c * yh_[jp1][i] + s * dky[i];
	}
	if(k == 0)
		return;
	r = pow(h_, (double)(-k));

	for(size_t i = 1; i <= n; i++)
		dky[i] *= r;

} /* end intdy   */

void LSODA::cfode(int meth_)
{
	int i, nq, nqm1, nqp1;
	double agamq, fnq, fnqm1, pc[13], pint, ragq, rqfac, rq1fac, tsign, xpin;
	/*
	   cfode is called by the integrator routine to set coefficients
	   needed there.  The coefficients for the current method, as
	   given by the value of meth_, are set for all orders and saved.
	   The maximum order assumed here is 12 if meth_ = 1 and 5 if meth_ = 2.
	   ( A smaller value of the maximum order is also allowed. )
	   cfode is called once at the beginning of the problem, and
	   is not called again unless and until meth_ is changed.

	   The elco array contains the basic method coefficients.
	   The coefficients el[i], 1 < i < nq+1, for the method of
	   order nq are stored in elco[nq][i].  They are given by a generating
	   polynomial, i.e.,

		  l(x) = el[1] + el[2]*x + ... + el[nq+1]*x^nq.

	   For the implicit Adams method, l(x) is given by

		  dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),   l(-1) = 0.

	   For the bdf methods, l(x) is given by

		  l(x) = (x+1)*(x+2)*...*(x+nq)/k,

	   where   k = factorial(nq)*(1+1/2+...+1/nq).

	   The tesco array contains test constants used for the
	   local error test and the selection of step size and/or order.
	   At order nq, tesco[nq][k] is used for the selection of step
	   size at order nq-1 if k = 1, at order nq if k = 2, and at order
	   nq+1 if k = 3.
	*/
	if(meth_ == 1) {
		elco[1][1]   = 1.;
		elco[1][2]   = 1.;
		tesco[1][1]  = 0.;
		tesco[1][2]  = 2.;
		tesco[2][1]  = 1.;
		tesco[12][3] = 0.;
		pc[1]        = 1.;
		rqfac        = 1.;
		for(nq = 2; nq <= 12; nq++) {
			/*
			   The pc array will contain the coefficients of the polynomial

				  p(x) = (x+1)*(x+2)*...*(x+nq-1).

			   Initially, p(x) = 1.
			*/
			rq1fac = rqfac;
			rqfac  = rqfac / (double)nq;
			nqm1   = nq - 1;
			fnqm1  = (double)nqm1;
			nqp1   = nq + 1;
			/*
			   Form coefficients of p(x)*(x+nq-1).
			*/
			pc[nq] = 0.;
			for(i = nq; i >= 2; i--)
				pc[i] = pc[i - 1] + fnqm1 * pc[i];
			pc[1] = fnqm1 * pc[1];
			/*
			   Compute integral, -1 to 0, of p(x) and x*p(x).
			*/
			pint  = pc[1];
			xpin  = pc[1] / 2.;
			tsign = 1.;
			for(i = 2; i <= nq; i++) {
				tsign = -tsign;
				pint += tsign * pc[i] / (double)i;
				xpin += tsign * pc[i] / (double)(i + 1);
			}
			/*
			   Store coefficients in elco and tesco.
			*/
			elco[nq][1] = pint * rq1fac;
			elco[nq][2] = 1.;
			for(i = 2; i <= nq; i++)
				elco[nq][i + 1] = rq1fac * pc[i] / (double)i;
			agamq        = rqfac * xpin;
			ragq         = 1. / agamq;
			tesco[nq][2] = ragq;
			if(nq < 12)
				tesco[nqp1][1] = ragq * rqfac / (double)nqp1;
			tesco[nqm1][3] = ragq;
		} /* end for   */
		return;
	} /* end if ( meth_ == 1 )   */

	/* meth_ = 2. */
	pc[1]  = 1.;
	rq1fac = 1.;

	/*
	   The pc array will contain the coefficients of the polynomial
		  p(x) = (x+1)*(x+2)*...*(x+nq).
	   Initially, p(x) = 1.
	*/
	for(nq = 1; nq <= 5; nq++) {
		fnq  = (double)nq;
		nqp1 = nq + 1;
		/*
		   Form coefficients of p(x)*(x+nq).
		*/
		pc[nqp1] = 0.;
		for(i = nq + 1; i >= 2; i--)
			pc[i] = pc[i - 1] + fnq * pc[i];
		pc[1] *= fnq;
		/*
		   Store coefficients in elco and tesco.
		*/
		for(i = 1; i <= nqp1; i++)
			elco[nq][i] = pc[i] / pc[2];
		elco[nq][2]  = 1.;
		tesco[nq][1] = rq1fac;
		tesco[nq][2] = ((double)nqp1) / elco[nq][1];
		tesco[nq][3] = ((double)(nq + 2)) / elco[nq][1];
		rq1fac /= fnq;
	}
	return;

} /* end cfode   */

void LSODA::scaleh(double *rh, double *pdh)
{
	double r;
	/*
	   If h_ is being changed, the h_ ratio rh is checked against rmax, hmin,
	   and hmxi, and the yh_ array is rescaled.  ialth is set to l = nq + 1
	   to prevent a change of h_ for that many steps, unless forced by a
	   convergence or error test failure.
	*/
	*rh = min(*rh, rmax);
	*rh = *rh / max(1., fabs(h_) * hmxi * *rh);
	/*
	   If meth_ = 1, also restrict the new step size by the stability region.
	   If this reduces h_, set irflag to 1 so that if there are roundoff
	   problems later, we can assume that is the cause of the trouble.
	*/
	if(meth_ == 1) {
		irflag = 0;
		*pdh   = max(fabs(h_) * pdlast, 0.000001);
		if((*rh * *pdh * 1.00001) >= sm1[nq]) {
			*rh    = sm1[nq] / *pdh;
			irflag = 1;
		}
	}
	r = 1.;
	for(size_t j = 2; j <= l; j++) {
		r *= *rh;
		for(size_t i = 1; i <= n; i++)
			yh_[j][i] *= r;
	}
	h_ *= *rh;
	rc *= *rh;
	ialth = l;

} /* end scaleh   */

/*
   This function routine computes the weighted max-norm
   of the vector of length n contained in the array v, with weights
   contained in the array w of length n.

   vmnorm = max( i = 1, ..., n ) fabs( v[i] ) * w[i].
*/
double LSODA::vmnorm(const size_t n, const vector<double> &v, const vector<double> &w)
{
	double vm = 0.;
	for(size_t i = 1; i <= n; i++)
		vm = max(vm, fabs(v[i]) * w[i]);
	return vm;
}

double LSODA::fnorm(int n, const vector<vector<double>> &a, const vector<double> &w)

/*
   This subroutine computes the norm of a full n by n matrix,
   stored in the array a, that is consistent with the weighted max-norm
   on vectors, with weights stored in the array w.

	  fnorm = max(i=1,...,n) ( w[i] * sum(j=1,...,n) fabs( a[i][j] ) / w[j] )
*/

{
	double an = 0, sum = 0;

	for(size_t i = 1; i <= (size_t)n; i++) {
		sum = 0.;
		for(size_t j = 1; j <= (size_t)n; j++)
			sum += fabs(a[i][j]) / w[j];
		an = max(an, sum * w[i]);
	}
	return an;
}


void LSODA::corfailure(double *told, double *rh, size_t *ncf, size_t *corflag)
{
	ncf++;
	rmax = 2.;
	tn_  = *told;
	for(size_t j = nq; j >= 1; j--)
		for(size_t i1 = j; i1 <= nq; i1++)
			for(size_t i = 1; i <= n; i++)
				yh_[i1][i] -= yh_[i1 + 1][i];

	if(fabs(h_) <= hmin * 1.00001 || *ncf == mxncf) {
		*corflag = 2;
		return;
	}
	*corflag = 1;
	*rh      = 0.25;
	ipup     = miter;
}

/*
   This routine manages the solution of the linear system arising from
   a chord iteration.  It is called if miter != 0.
   If miter is 2, it calls dgesl to accomplish this.
   If miter is 5, it calls dgbsl.

   y = the right-hand side vector on input, and the solution vector
	   on output.
*/
void LSODA::solsy(vector<double> &y)
{
	iersl = 0;
	if(miter != 2) {
		printf("solsy -- miter != 2\n");
		return;
	}
	if(miter == 2)
		dgesl(wm_, n, ipvt, y, 0);
	return;
}

void LSODA::methodswitch(double dsm, double pnorm, double *pdh, double *rh)
{
	int lm1, lm1p1, lm2, lm2p1, nqm1, nqm2;
	double rh1, rh2, rh1it, exm2, dm2, exm1, dm1, alpha, exsm;

	/*
	   We are current using an Adams method.  Consider switching to bdf.
	   If the current order is greater than 5, assume the problem is
	   not stiff, and skip this section.
	   If the Lipschitz constant and error estimate are not polluted
	   by roundoff, perform the usual test.
	   Otherwise, switch to the bdf methods if the last step was
	   restricted to insure stability ( irflag = 1 ), and stay with Adams
	   method if not.  When switching to bdf with polluted error estimates,
	   in the absence of other information, double the step size.

	   When the estimates are ok, we make the usual test by computing
	   the step size we could have (ideally) used on this step,
	   with the current (Adams) method, and also that for the bdf.
	   If nq > mxords, we consider changing to order mxords on switching.
	   Compare the two step sizes to decide whether to switch.
	   The step size advantage must be at least ratio = 5 to switch.
	*/
	if(meth_ == 1) {
		if(nq > 5)
			return;
		if(dsm <= (100. * pnorm * ETA) || pdest == 0.) {
			if(irflag == 0)
				return;
			rh2  = 2.;
			nqm2 = min(nq, mxords);
		}
		else {
			exsm  = 1. / (double)l;
			rh1   = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);
			rh1it = 2. * rh1;
			*pdh  = pdlast * fabs(h_);
			if((*pdh * rh1) > 0.00001)
				rh1it = sm1[nq] / *pdh;
			rh1 = min(rh1, rh1it);
			if(nq > mxords) {
				nqm2  = mxords;
				lm2   = mxords + 1;
				exm2  = 1. / (double)lm2;
				lm2p1 = lm2 + 1;
				dm2   = vmnorm(n, yh_[lm2p1], ewt) / cm2[mxords];
				rh2   = 1. / (1.2 * pow(dm2, exm2) + 0.0000012);
			}
			else {
				dm2  = dsm * (cm1[nq] / cm2[nq]);
				rh2  = 1. / (1.2 * pow(dm2, exsm) + 0.0000012);
				nqm2 = nq;
			}
			if(rh2 < ratio * rh1)
				return;
		}
		/*
		   The method switch test passed.  Reset relevant quantities for bdf.
		*/
		*rh    = rh2;
		icount = 20;
		meth_  = 2;
		miter  = jtyp;
		pdlast = 0.;
		nq     = nqm2;
		l      = nq + 1;
		return;
	} /* end if ( meth_ == 1 )   */

	/*
	   We are currently using a bdf method, considering switching to Adams.
	   Compute the step size we could have (ideally) used on this step,
	   with the current (bdf) method, and also that for the Adams.
	   If nq > mxordn, we consider changing to order mxordn on switching.
	   Compare the two step sizes to decide whether to switch.
	   The step size advantage must be at least 5/ratio = 1 to switch.
	   If the step size for Adams would be so small as to cause
	   roundoff pollution, we stay with bdf.
	*/
	exsm = 1. / (double)l;
	if(mxordn < nq) {
		nqm1  = mxordn;
		lm1   = mxordn + 1;
		exm1  = 1. / (double)lm1;
		lm1p1 = lm1 + 1;
		dm1   = vmnorm(n, yh_[lm1p1], ewt) / cm1[mxordn];
		rh1   = 1. / (1.2 * pow(dm1, exm1) + 0.0000012);
	}
	else {
		dm1  = dsm * (cm2[nq] / cm1[nq]);
		rh1  = 1. / (1.2 * pow(dm1, exsm) + 0.0000012);
		nqm1 = nq;
		exm1 = exsm;
	}
	rh1it = 2. * rh1;
	*pdh  = pdnorm * fabs(h_);
	if((*pdh * rh1) > 0.00001)
		rh1it = sm1[nqm1] / *pdh;
	rh1 = min(rh1, rh1it);
	rh2 = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);
	if((rh1 * ratio) < (5. * rh2))
		return;
	alpha = max(0.001, rh1);
	dm1 *= pow(alpha, exm1);
	if(dm1 <= 1000. * ETA * pnorm)
		return;
	/*
	   The switch test passed.  Reset relevant quantities for Adams.
	*/
	*rh    = rh1;
	icount = 20;
	meth_  = 1;
	miter  = 0;
	pdlast = 0.;
	nq     = nqm1;
	l      = nq + 1;
} /* end methodswitch   */

/*
   This routine returns from stoda to lsoda.  Hence freevectors() is
   not executed.
*/
void LSODA::endstoda()
{
	double r = 1. / tesco[nqu][2];
	for(size_t i = 1; i <= n; i++)
		acor[i] *= r;
	hold   = h_;
	jstart = 1;
}

/*
   Regardless of the success or failure of the step, factors
   rhdn, rhsm, and rhup are computed, by which h_ could be multiplied
   at order nq - 1, order nq, or order nq + 1, respectively.
   In the case of a failure, rhup = 0. to avoid an order increase.
   The largest of these is determined and the new order chosen
   accordingly.  If the order is to be increased, we compute one
   additional scaled derivative.

   orderflag = 0  : no change in h_ or nq,
			   1  : change in h_ but not nq,
			   2  : change in both h_ and nq.
*/
void LSODA::orderswitch(
	double *rhup, double dsm, double *pdh, double *rh, size_t *orderflag)
{
	size_t newq = 0;
	double exsm, rhdn, rhsm, ddn, exdn, r;

	*orderflag = 0;

	exsm = 1. / (double)l;
	rhsm = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);

	rhdn = 0.;
	if(nq != 1) {
		ddn  = vmnorm(n, yh_[l], ewt) / tesco[nq][1];
		exdn = 1. / (double)nq;
		rhdn = 1. / (1.3 * pow(ddn, exdn) + 0.0000013);
	}
	/*
	   If meth_ = 1, limit rh accordinfg to the stability region also.
	*/
	if(meth_ == 1) {
		*pdh = max(fabs(h_) * pdlast, 0.000001);
		if(l < lmax)
			*rhup = min(*rhup, sm1[l] / *pdh);
		rhsm = min(rhsm, sm1[nq] / *pdh);
		if(nq > 1)
			rhdn = min(rhdn, sm1[nq - 1] / *pdh);
		pdest = 0.;
	}
	if(rhsm >= *rhup) {
		if(rhsm >= rhdn) {
			newq = nq;
			*rh  = rhsm;
		}
		else {
			newq = nq - 1;
			*rh  = rhdn;
			if(kflag < 0 && *rh > 1.)
				*rh = 1.;
		}
	}
	else {
		if(*rhup <= rhdn) {
			newq = nq - 1;
			*rh  = rhdn;
			if(kflag < 0 && *rh > 1.)
				*rh = 1.;
		}
		else {
			*rh = *rhup;
			if(*rh >= 1.1) {
				r  = el[l] / (double)l;
				nq = l;
				l  = nq + 1;
				for(size_t i = 1; i <= n; i++)
					yh_[l][i] = acor[i] * r;

				*orderflag = 2;
				return;
			}
			else {
				ialth = 3;
				return;
			}
		}
	}
	/*
	   If meth_ = 1 and h_ is restricted by stability, bypass 10 percent test.
	*/
	if(1 == meth_) {
		if((*rh * *pdh * 1.00001) < sm1[newq])
			if(kflag == 0 && *rh < 1.1) {
				ialth = 3;
				return;
			}
	}
	else {
		if(kflag == 0 && *rh < 1.1) {
			ialth = 3;
			return;
		}
	}
	if(kflag <= -2)
		*rh = min(*rh, 0.2);
	/*
	   If there is a change of order, reset nq, l, and the coefficients.
	   In any case h_ is reset according to rh and the yh_ array is rescaled.
	   Then exit or redo the step.
	*/
	if(newq == nq) {
		*orderflag = 1;
		return;
	}
	nq         = newq;
	l          = nq + 1;
	*orderflag = 2;

} /* end orderswitch   */

void LSODA::resetcoeff()
/*
   The el vector and related constants are reset
   whenever the order nq is changed, or at the start of the problem.
*/
{
	array<double, 14> ep1;

	ep1 = elco[nq];
	for(size_t i = 1; i <= l; i++)
		el[i] = ep1[i];
	rc    = rc * el[1] / el0;
	el0   = el[1];
	conit = 0.5 / (double)(nq + 2);
}

void LSODA::_freevectors(void)
{
	// Does nothing. USE c++ memory mechanism here.
}


void LSODA::resize_system(int nyh, int lenyh){
	// yh_ and wm_ need a clear() because otherwise, only additional elements will have the new length nyh+1
	yh_.clear(); yh_.resize(lenyh + 1, std::vector<double>(nyh + 1, 0.0));
	wm_.clear(); wm_.resize(nyh + 1, std::vector<double>(nyh + 1, 0.0));
	ewt.resize(1 + nyh, 0);
	savf.resize(1 + nyh, 0);
	acor.resize(nyh + 1, 0.0);
	ipvt.resize(nyh + 1, 0.0);
}


int LSODA::get_fncalls(){
	return nfe;
}


int LSODA::get_istate() const{
	return iState;
}

void LSODA::set_istate(int n){
	iState = n;
}
