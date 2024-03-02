#include <iostream>
#include <cassert>

	
template <class Functor>
void LSODA::prja(
	const size_t neq, std::vector<double> &y, Functor &derivs, void *_data)
{
	(void)neq;

	size_t i = 0, ier = 0, j = 0;
	double fac = 0.0, hl0 = 0.0, r = 0.0, r0 = 0.0, yj = 0.0;
	/*
	   prja is called by stoda to compute and process the matrix
	   P = I - h_ * el[1] * J, where J is an approximation to the Jacobian.
	   Here J is computed by finite differencing.
	   J, scaled by -h_ * el[1], is stored in wm_.  Then the norm of J ( the
	   matrix norm consistent with the weighted max-norm on vectors given
	   by vmnorm ) is computed, and J is overwritten by P.  P is then
	   subjected to LU decomposition in preparation for later solution
	   of linear systems with p as coefficient matrix.  This is done
	   by dgefa if miter = 2, and by dgbfa if miter = 5.
	*/
	nje++;
	ierpj = 0;
	jcur  = 1;
	hl0   = h_ * el0;
	/*
	   If miter = 2, make n calls to f to approximate J.
	*/
	if(miter != 2) {
		fprintf(stderr, "[prja] miter != 2\n");
		return;
	}
	if(miter == 2) {
		fac = vmnorm(n, savf, ewt);
		r0  = 1000. * fabs(h_) * ETA * ((double)n) * fac;
		if(r0 == 0.)
			r0 = 1.;
		for(j = 1; j <= n; j++) {
			yj = y[j];
			r  = std::max(sqrteta * fabs(yj), r0 / ewt[j]);
			y[j] += r;
			fac = -hl0 / r;
			derivs(tn_, ++y.begin(), ++acor.begin(), _data);
			for(i = 1; i <= n; i++)
				wm_[i][j] = (acor[i] - savf[i]) * fac;
			y[j] = yj;
		}
		nfe += n;
		/*
		   Compute norm of Jacobian.
		*/
		pdnorm = fnorm(n, wm_, ewt) / fabs(hl0);
		/*
		   Add identity matrix.
		*/
		for(i = 1; i <= n; i++)
			wm_[i][i] += 1.;
		/*
		   Do LU decomposition on P.
		*/
		dgefa(wm_, n, ipvt, &ier);
		if(ier != 0)
			ierpj = 1;
		return;
	}
} /* end prja   */



/*
 *corflag = 0 : corrector converged,
1 : step size to be reduced, redo prediction,
2 : corrector cannot converge, failure flag.
*/
template <class Functor>
void LSODA::correction(const size_t neq, std::vector<double> &y, Functor &derivs,
	size_t *corflag, double pnorm, double *del, double *delp, double *told, size_t *ncf,
	double *rh, size_t *m, void *_data)
{
	double rm = 0.0, rate = 0.0, dcon = 0.0;

	/*
	   Up to maxcor corrector iterations are taken.  A convergence test is
	   made on the r.m.s. norm of each correction, weighted by the error
	   weight vector ewt.  The sum of the corrections is accumulated in the
	   vector acor[i].  The yh_ array is not altered in the corrector loop.
	*/

	*m       = 0;
	*corflag = 0;
	*del     = 0.;

	for(size_t i = 1; i <= n; i++)
		y[i] = yh_[1][i];

	derivs(tn_, ++y.begin(), ++savf.begin(), _data);

	nfe++;
	/*
	   If indicated, the matrix P = I - h_ * el[1] * J is reevaluated and
	   preprocessed before starting the corrector iteration.  ipup is set
	   to 0 as an indicator that this has been done.
	*/
	while(1) {
		if(*m == 0) {
			if(ipup > 0) {
				prja(neq, y, derivs, _data);
				ipup  = 0;
				rc    = 1.;
				nslp  = nst;
				crate = 0.7;
				if(ierpj != 0) {
					corfailure(told, rh, ncf, corflag);
					return;
				}
			}
			for(size_t i = 1; i <= n; i++)
				acor[i] = 0.;
		} /* end if ( *m == 0 )   */
		if(miter == 0) {
			/*
			   In case of functional iteration, update y directly from
			   the result of the last function evaluation.
			*/
			for(size_t i = 1; i <= n; i++) {
				savf[i] = h_ * savf[i] - yh_[2][i];
				y[i]    = savf[i] - acor[i];
			}
			*del = vmnorm(n, y, ewt);
			for(size_t i = 1; i <= n; i++) {
				y[i]    = yh_[1][i] + el[1] * savf[i];
				acor[i] = savf[i];
			}
		}
		/* end functional iteration   */
		/*
		   In the case of the chord method, compute the corrector error,
		   and solve the linear system with that as right-hand side and
		   P as coefficient matrix.
		 */
		else {
			for(size_t i = 1; i <= n; i++)
				y[i] = h_ * savf[i] - (yh_[2][i] + acor[i]);

			solsy(y);
			*del = vmnorm(n, y, ewt);

			for(size_t i = 1; i <= n; i++) {
				acor[i] += y[i];
				y[i] = yh_[1][i] + el[1] * acor[i];
			}
		} /* end chord method   */
		/*
		   Test for convergence.  If *m > 0, an estimate of the convergence
		   rate constant is stored in crate, and this is used in the test.

		   We first check for a change of iterates that is the size of
		   roundoff error.  If this occurs, the iteration has converged, and a
		   new rate estimate is not formed.
		   In all other cases, force at least two iterations to estimate a
		   local Lipschitz constant estimate for Adams method.
		   On convergence, form pdest = local maximum Lipschitz constant
		   estimate.  pdlast is the most recent nonzero estimate.
		*/
		if(*del <= 100. * pnorm * ETA)
			break;
		if(*m != 0 || meth_ != 1) {
			if(*m != 0) {
				rm = 1024.0;
				if(*del <= (1024. * *delp))
					rm = *del / *delp;
				rate  = std::max(rate, rm);
				crate = std::max(0.2 * crate, rm);
			}
			dcon = *del * std::min(1., 1.5 * crate) / (tesco[nq][2] * conit);
			if(dcon <= 1.) {
				pdest = std::max(pdest, rate / fabs(h_ * el[1]));
				if(pdest != 0.)
					pdlast = pdest;
				break;
			}
		}
		/*
		   The corrector iteration failed to converge.
		   If miter != 0 and the Jacobian is out of date, prja is called for
		   the next try.   Otherwise the yh_ array is retracted to its values
		   before prediction, and h_ is reduced, if possible.  If h_ cannot be
		   reduced or mxncf failures have occured, exit with corflag = 2.
		*/
		(*m)++;
		if(*m == maxcor || (*m >= 2 && *del > 2. * *delp)) {
			if(miter == 0 || jcur == 1) {
				corfailure(told, rh, ncf, corflag);
				return;
			}
			ipup = miter;
			/*
			   Restart corrector if Jacobian is recomputed.
			*/
			*m   = 0;
			rate = 0.;
			*del = 0.;
			for(size_t i = 1; i <= n; i++)
				y[i] = yh_[1][i];

			derivs(tn_, ++y.begin(), ++savf.begin(), _data);

			nfe++;
		}
		/*
		   Iterate corrector.
		*/
		else {
			*delp = *del;
			derivs(tn_, ++y.begin(), ++savf.begin(), _data);
			nfe++;
		}
	} /* end while   */
} /* end correction   */


template <class Functor>
void LSODA::stoda(
	const size_t neq, std::vector<double> &y, Functor &derivs, void *_data)
{
	assert(neq + 1 == y.size());

	size_t corflag = 0, orderflag = 0;
	size_t i = 0, i1 = 0, j = 0, m = 0, ncf = 0;
	double del = 0.0, delp = 0.0, dsm = 0.0, dup = 0.0, exup = 0.0, r = 0.0, rh = 0.0,
		   rhup = 0.0, told = 0.0;
	double pdh = 0.0, pnorm = 0.0;

	/*
	   stoda performs one step of the integration of an initial value
	   problem for a system of ordinary differential equations.
	   Note.. stoda is independent of the value of the iteration method
	   indicator miter, when this is != 0, and hence is independent
	   of the type of chord method used, or the Jacobian structure.
	   Communication with stoda is done with the following variables:

	   jstart = an integer used for input only, with the following
				values and meanings:

				   0  perform the first step,
				 > 0  take a new step continuing from the last,
				  -1  take the next step with a new value of h_,
					  n, meth_, miter, and/or matrix parameters.
				  -2  take the next step with a new value of h_,
					  but with other inputs unchanged.

	   kflag = a completion code with the following meanings:

				 0  the step was successful,
				-1  the requested error could not be achieved,
				-2  corrector convergence could not be achieved,
				-3  fatal error in prja or solsy.

	   miter = corrector iteration method:

				 0  functional iteration,
				>0  a chord method corresponding to jacobian type jt.

	*/
	kflag = 0;
	told  = tn_;
	ncf   = 0;
	ierpj = 0;
	iersl = 0;
	jcur  = 0;
	delp  = 0.;

	/*
	   On the first call, the order is set to 1, and other variables are
	   initialized.  rmax is the maximum ratio by which h_ can be increased
	   in a single step.  It is initially 1.e4 to compensate for the small
	   initial h_, but then is normally equal to 10.  If a filure occurs
	   (in corrector convergence or error test), rmax is set at 2 for
	   the next increase.
	   cfode is called to get the needed coefficients for both methods.
	*/
	if(jstart == 0) {
		lmax  = maxord + 1;
		nq    = 1;
		l     = 2;
		ialth = 2;
		rmax  = 10000.;
		rc    = 0.;
		el0   = 1.;
		crate = 0.7;
		hold  = h_;
		nslp  = 0;
		ipup  = miter;
		/*
		   Initialize switching parameters.  meth_ = 1 is assumed initially.
		*/
		icount = 20;
		irflag = 0;
		pdest  = 0.;
		pdlast = 0.;
		ratio  = 5.;
		cfode(2);
		for(i = 1; i <= 5; i++)
			cm2[i] = tesco[i][2] * elco[i][i + 1];
		cfode(1);
		for(i = 1; i <= 12; i++)
			cm1[i] = tesco[i][2] * elco[i][i + 1];
		resetcoeff();
	} /* end if ( jstart == 0 )   */
	/*
	   The following block handles preliminaries needed when jstart = -1.
	   ipup is set to miter to force a matrix update.
	   If an order increase is about to be considered ( ialth = 1 ),
	   ialth is reset to 2 to postpone consideration one more step.
	   If the caller has changed meth_, cfode is called to reset
	   the coefficients of the method.
	   If h_ is to be changed, yh_ must be rescaled.
	   If h_ or meth_ is being changed, ialth is reset to l = nq + 1
	   to prevent further changes in h_ for that many steps.
	*/
	if(jstart == -1) {
		ipup = miter;
		lmax = maxord + 1;
		if(ialth == 1)
			ialth = 2;
		if(meth_ != mused) {
			cfode(meth_);
			ialth = l;
			resetcoeff();
		}
		if(h_ != hold) {
			rh = h_ / hold;
			h_ = hold;
			scaleh(&rh, &pdh);
		}
	} /* if ( jstart == -1 )   */
	if(jstart == -2) {
		if(h_ != hold) {
			rh = h_ / hold;
			h_ = hold;
			scaleh(&rh, &pdh);
		}
	} /* if ( jstart == -2 )   */

	/*
	   Prediction.
	   This section computes the predicted values by effectively
	   multiplying the yh_ array by the pascal triangle matrix.
	   rc is the ratio of new to old values of the coefficient h_ * el[1].
	   When rc differs from 1 by more than ccmax, ipup is set to miter
	   to force pjac to be called, if a jacobian is involved.
	   In any case, prja is called at least every msbp steps.
	*/
	while(1) {
		while(1) {
			if(fabs(rc - 1.) > ccmax)
				ipup = miter;
			if(nst >= nslp + msbp)
				ipup = miter;
			tn_ += h_;
			for(size_t j = nq; j >= 1; j--)
				for(size_t i1 = j; i1 <= nq; i1++)
					for(i = 1; i <= n; i++)
						yh_[i1][i] += yh_[i1 + 1][i];

			pnorm = vmnorm(n, yh_[1], ewt);
			correction(
				neq, y, derivs, &corflag, pnorm, &del, &delp, &told, &ncf, &rh, &m, _data);
			if(corflag == 0)
				break;
			if(corflag == 1) {
				rh = std::max(rh, hmin / fabs(h_));
				scaleh(&rh, &pdh);
				continue;
			}
			if(corflag == 2) {
				kflag  = -2;
				hold   = h_;
				jstart = 1;
				return;
			}
		} /* end inner while ( corrector loop )   */

		/*
		   The corrector has converged.  jcur is set to 0
		   to signal that the Jacobian involved may need updating later.
		   The local error test is done now.
		*/
		jcur = 0;
		if(m == 0)
			dsm = del / tesco[nq][2];
		if(m > 0)
			dsm = vmnorm(n, acor, ewt) / tesco[nq][2];

		if(dsm <= 1.) {
			/*
			   After a successful step, update the yh_ array.
			   Decrease icount by 1, and if it is -1, consider switching methods.
			   If a method switch is made, reset various parameters,
			   rescale the yh_ array, and exit.  If there is no switch,
			   consider changing h_ if ialth = 1.  Otherwise decrease ialth by 1.
			   If ialth is then 1 and nq < maxord, then acor is saved for
			   use in a possible order increase on the next step.
			   If a change in h_ is considered, an increase or decrease in order
			   by one is considered also.  A change in h_ is made only if it is by
			   a factor of at least 1.1.  If not, ialth is set to 3 to prevent
			   testing for that many steps.
			*/
			kflag = 0;
			nst++;
			hu    = h_;
			nqu   = nq;
			mused = meth_;
			for(size_t j = 1; j <= l; j++) {
				r = el[j];
				for(i = 1; i <= n; i++)
					yh_[j][i] += r * acor[i];
			}
			icount--;
			if(icount < 0) {
				methodswitch(dsm, pnorm, &pdh, &rh);
				if(meth_ != mused) {
					rh = std::max(rh, hmin / fabs(h_));
					scaleh(&rh, &pdh);
					rmax = 10.;
					endstoda();
					break;
				}
			}
			/*
			   No method switch is being made.  Do the usual step/order selection.
			*/
			ialth--;
			if(ialth == 0) {
				rhup = 0.;
				if(l != lmax) {
					for(i = 1; i <= n; i++)
						savf[i] = acor[i] - yh_[lmax][i];
					dup  = vmnorm(n, savf, ewt) / tesco[nq][3];
					exup = 1. / (double)(l + 1);
					rhup = 1. / (1.4 * pow(dup, exup) + 0.0000014);
				}

				orderswitch(&rhup, dsm, &pdh, &rh, &orderflag);

				/*
				   No change in h_ or nq.
				*/
				if(orderflag == 0) {
					endstoda();
					break;
				}
				/*
				   h_ is changed, but not nq.
				*/
				if(orderflag == 1) {
					rh = std::max(rh, hmin / fabs(h_));
					scaleh(&rh, &pdh);
					rmax = 10.;
					endstoda();
					break;
				}
				/*
				   both nq and h_ are changed.
				*/
				if(orderflag == 2) {
					resetcoeff();
					rh = std::max(rh, hmin / fabs(h_));
					scaleh(&rh, &pdh);
					rmax = 10.;
					endstoda();
					break;
				}
			} /* end if ( ialth == 0 )   */
			if(ialth > 1 || l == lmax) {
				endstoda();
				break;
			}

			for(size_t i = 1; i <= n; i++)
				yh_[lmax][i] = acor[i];

			endstoda();
			break;
		}
		/* end if ( dsm <= 1. )   */
		/*
		   The error test failed.  kflag keeps track of multiple failures.
		   Restore tn_ and the yh_ array to their previous values, and prepare
		   to try the step again.  Compute the optimum step size for this or
		   one lower.  After 2 or more failures, h_ is forced to decrease
		   by a factor of 0.2 or less.
		 */
		else {
			kflag--;
			tn_ = told;
			for(j = nq; j >= 1; j--) {
				for(i1 = j; i1 <= nq; i1++)
					for(i = 1; i <= n; i++)
						yh_[i1][i] -= yh_[i1 + 1][i];
			}
			rmax = 2.;
			if(fabs(h_) <= hmin * 1.00001) {
				kflag  = -1;
				hold   = h_;
				jstart = 1;
				break;
			}
			if(kflag > -3) {
				rhup = 0.;
				orderswitch(&rhup, dsm, &pdh, &rh, &orderflag);
				if(orderflag == 1 || orderflag == 0) {
					if(orderflag == 0)
						rh = std::min(rh, 0.2);
					rh = std::max(rh, hmin / fabs(h_));
					scaleh(&rh, &pdh);
				}
				if(orderflag == 2) {
					resetcoeff();
					rh = std::max(rh, hmin / fabs(h_));
					scaleh(&rh, &pdh);
				}
				continue;
			}
			/* if ( kflag > -3 )   */
			/*
			   Control reaches this section if 3 or more failures have occurred.
			   If 10 failures have occurred, exit with kflag = -1.
			   It is assumed that the derivatives that have accumulated in the
			   yh_ array have errors of the wrong order.  Hence the first
			   derivative is recomputed, and the order is set to 1.  Then
			   h_ is reduced by a factor of 10, and the step is retried,
			   until it succeeds or h_ reaches hmin.
			 */
			else {
				if(kflag == -10) {
					kflag  = -1;
					hold   = h_;
					jstart = 1;
					break;
				}
				else {
					rh = 0.1;
					rh = std::max(hmin / fabs(h_), rh);
					h_ *= rh;
					for(i = 1; i <= n; i++)
						y[i] = yh_[1][i];
					derivs(tn_, ++y.begin(), ++savf.begin(), _data);
					nfe++;
					for(i = 1; i <= n; i++)
						yh_[2][i] = h_ * savf[i];
					ipup  = miter;
					ialth = 5;
					if(nq == 1)
						continue;
					nq = 1;
					l  = 2;
					resetcoeff();
					continue;
				}
			} /* end else -- kflag <= -3 */
		}     /* end error failure handling   */
	}         /* end outer while   */

} /* end stoda   */


/*
   The following block handles all successful returns from lsoda.
   If itask != 1, y is loaded from yh_ and t is set accordingly.
   *Istate is set to 2, the illegal input counter is zeroed, and the
   optional outputs are loaded into the work arrays before returning.
*/

template<class AfterStep>
void LSODA::successreturn(AfterStep &after_step,
	std::vector<double> &y, double *t, int itask, int ihit, double tcrit, int *istate)
{
	for(size_t i = 1; i <= n; i++)
		y[i] = yh_[1][i];
	*t = tn_;
	if(itask == 4 || itask == 5)
		if(ihit)
			*t = tcrit;
	*istate = 2;
	illin   = 0;
	
	after_step(*t, ++y.begin());
}



/*
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c 2.  linda r. petzold, automatic selection of methods for solving
c     stiff and nonstiff systems of ordinary differential equations,
c     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
c-----------------------------------------------------------------------
*/
template <class Functor, class AfterStep>
void LSODA::lsoda(Functor &derivs, AfterStep &after_step, const size_t neq, std::vector<double> &y, double *t,
	double tout, int itask, int *istate, int iopt, int jt, std::array<int, 7> &iworks,
	std::array<double, 4> &rworks, void *_data)
{
	assert(tout > *t);

	int mxstp0 = 500, mxhnl0 = 10;

	int iflag = 0, lenyh = 0, ihit = 0;

	double atoli = 0, ayi = 0, big = 0, h0 = 0, hmax = 0, hmx = 0, rh = 0, rtoli = 0,
		   tcrit = 0, tdist = 0, tnext = 0, tol = 0, tolsf = 0, tp = 0, size = 0, sum = 0,
		   w0 = 0;

	/*
	   Block a.
	   This code block is executed on every call.
	   It tests *istate and itask for legality and branches appropriately.
	   If *istate > 1 but the flag init shows that initialization has not
	   yet been done, an error return occurs.
	   If *istate = 1 and tout = t, return immediately.
	*/

	if(*istate < 1 || *istate > 3) {
		// fprintf(stderr, "[lsoda] illegal istate = %d\n", *istate);
		std::cerr << "[lsoda] illegal istate = " << *istate << std::endl;
		terminate(istate);
		return;
	}
	if(itask < 1 || itask > 5) {
		fprintf(stderr, "[lsoda] illegal itask = %d\n", itask);
		terminate(istate);
		return;
	}
	if(init == 0 && (*istate == 2 || *istate == 3)) {
		fprintf(stderr, "[lsoda] istate > 1 but lsoda not initialized\n");
		terminate(istate);
		return;
	}

	/*
	   Block b.
	   The next code block is executed for the initial call ( *istate = 1 ),
	   or for a continuation call with parameter changes ( *istate = 3 ).
	   It contains checking of all inputs and various initializations.

	   First check legality of the non-optional inputs neq, itol, iopt,
	   jt, ml, and mu.
	*/

	if(*istate == 1 || *istate == 3) {
		ntrep = 0;
		if(neq <= 0) {
			std::cerr << "[lsoda] neq = " << neq << " is less than 1." << std::endl;
			terminate(istate);
			return;
		}
		if(*istate == 3 && neq > n) {
			std::cerr << "[lsoda] istate = 3 and neq increased" << std::endl;
			terminate(istate);
			return;
		}
		n = neq;
		if(itol_ < 1 || itol_ > 4) {
			std::cerr << "[lsoda] itol = " << itol_ << " illegal" << std::endl;
			terminate(istate);
			return;
		}
		if(iopt < 0 || iopt > 1) {
			std::cerr << "[lsoda] iopt = " << iopt << " illegal" << std::endl;
			terminate(istate);
			return;
		}
		if(jt == 3 || jt < 1 || jt > 5) {
			std::cerr << "[lsoda] jt = " << jt << " illegal" << std::endl;
			terminate(istate);
			return;
		}
		jtyp = jt;
		if(jt > 2) {
			ml = iworks[0];
			mu = iworks[1];
			if(ml >= n) {
				std::cerr << "[lsoda] ml = " << ml << " not between 1 and neq" << std::endl;
				terminate(istate);
				return;
			}
			if(mu >= n) {
				std::cerr << "[lsoda] mu = " << mu << " not between 1 and neq" << std::endl;
				terminate(istate);
				return;
			}
		}

		/* Next process and check the optional inpus.   */
		/* Default options.   */
		if(iopt == 0) {
			ixpr   = 0;
			mxstep = mxstp0;
			mxhnil = mxhnl0;
			hmxi   = 0.;
			hmin   = 0.;
			if(*istate == 1) {
				h0     = 0.;
				mxordn = mord[0];
				mxords = mord[1];
			}
		}
		/* end if ( iopt == 0 )   */
		/* Optional inputs.   */
		else /* if ( iopt = 1 )  */
		{
			ixpr = iworks[2];
			if(ixpr > 1) {
				std::cerr << "[lsoda] ixpr = " << ixpr << " is illegal" << std::endl;
				terminate(istate);
				return;
			}

			mxstep = iworks[3];
			if(mxstep == 0)
				mxstep = mxstp0;
			mxhnil = iworks[4];

			if(*istate == 1) {
				h0     = rworks[1];
				mxordn = iworks[5];

				if(mxordn == 0)
					mxordn = 100;

				mxordn = std::min(mxordn, mord[0]);
				mxords = iworks[6];

				// if mxords is not given use 100.
				if(mxords == 0)
					mxords = 100;

				mxords = std::min(mxords, mord[1]);

				if((tout - *t) * h0 < 0.) {
					std::cerr << "[lsoda] tout = " << tout << " behind t = " << *t
						 << ". integration direction is given by " << h0 << std::endl;
					terminate(istate);
					return;
				}
			} /* end if ( *istate == 1 )  */
			hmax = rworks[2];
			if(hmax < 0.) {
				std::cerr << "[lsoda] hmax < 0." << std::endl;
				terminate(istate);
				return;
			}
			hmxi = 0.;
			if(hmax > 0)
				hmxi = 1. / hmax;

			hmin = rworks[3];
			if(hmin < 0.) {
				std::cerr << "[lsoda] hmin < 0." << std::endl;
				terminate(istate);
				return;
			}
		} /* end else   */ /* end iopt = 1   */
	}                      /* end if ( *istate == 1 || *istate == 3 )   */
	/*
	   If *istate = 1, meth_ is initialized to 1.

	   Also allocate memory for yh_, wm_, ewt, savf, acor, ipvt.
	*/
	if(*istate == 1) {
		/*
		   If memory were not freed, *istate = 3 need not reallocate memory.
		   Hence this section is not executed by *istate = 3.
		*/
		sqrteta = sqrt(ETA);
		meth_   = 1;

		nyh   = n;
		lenyh = 1 + std::max(mxordn, mxords);
		resize_system(nyh, lenyh);

		//// yh_ and wm_ need a clear() because otherwise, only additional elements will have the new length nyh
		//yh_.clear(); yh_.resize(lenyh + 1, std::vector<double>(nyh + 1, 0.0));
		//wm_.clear(); wm_.resize(nyh + 1, std::vector<double>(nyh + 1, 0.0));
		//ewt.resize(1 + nyh, 0);
		//savf.resize(1 + nyh, 0);
		//acor.resize(nyh + 1, 0.0);
		//ipvt.resize(nyh + 1, 0.0);
	}
	/*
	   Check rtol and atol for legality.
	*/
	if(*istate == 1 || *istate == 3) {
		rtoli = rtol_[1];
		atoli = atol_[1];
		for(size_t i = 1; i <= n; i++) {
			if(itol_ >= 3)
				rtoli = rtol_[i];
			if(itol_ == 2 || itol_ == 4)
				atoli = atol_[i];
			if(rtoli < 0.) {
				fprintf(stderr, "[lsoda] rtol = %g is less than 0.\n", rtoli);
				terminate(istate);
				return;
			}
			if(atoli < 0.) {
				fprintf(stderr, "[lsoda] atol = %g is less than 0.\n", atoli);
				terminate(istate);
				return;
			}
		} /* end for   */
	}     /* end if ( *istate == 1 || *istate == 3 )   */

	/* If *istate = 3, set flag to signal parameter changes to stoda. */
	if(*istate == 3) {
		jstart = -1;
	}
	/*
	   Block c.
	   The next block is for the initial call only ( *istate = 1 ).
	   It contains all remaining initializations, the initial call to f,
	   and the calculation of the initial step size.
	   The error weights in ewt are inverted after being loaded.
	*/
	if(*istate == 1) {
		tn_    = *t;
		tsw    = *t;
		maxord = mxordn;
		if(itask == 4 || itask == 5) {
			tcrit = rworks[0];
			if((tcrit - tout) * (tout - *t) < 0.) {
				fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tout\n");
				terminate(istate);
				return;
			}
			if(h0 != 0. && (*t + h0 - tcrit) * h0 > 0.)
				h0 = tcrit - *t;
		}

		jstart = 0;
		nhnil  = 0;
		nst    = 0;
		nje    = 0;
		nslast = 0;
		hu     = 0.;
		nqu    = 0;
		mused  = 0;
		miter  = 0;
		ccmax  = 0.3;
		maxcor = 3;
		msbp   = 20;
		mxncf  = 10;

		/* Initial call to f.  */
		assert((int)yh_.size() == lenyh + 1);
		assert(yh_[0].size() == nyh + 1);

		derivs(*t, ++y.begin(), ++yh_[2].begin(), _data);
		nfe = 1;

		/* Load the initial value vector in yh_.  */
		for(size_t i = 1; i <= n; i++)
			yh_[1][i] = y[i];

		/* Load and invert the ewt array.  ( h_ is temporarily set to 1. ) */
		nq = 1;
		h_ = 1.;
		ewset(y);
		for(size_t i = 1; i <= n; i++) {
			if(ewt[i] <= 0.) {
				std::cerr << "[lsoda] ewt[" << i << "] = " << ewt[i] << " <= 0.\n" << std::endl;
				terminate2(y, t);
				return;
			}
			ewt[i] = 1. / ewt[i];
		}

		/*
		   The coding below computes the step size, h0, to be attempted on the
		   first step, unless the user has supplied a value for this.
		   First check that tout - *t differs significantly from zero.
		   A scalar tolerance quantity tol is computed, as max(rtol[i])
		   if this is positive, or max(atol[i]/fabs(y[i])) otherwise, adjusted
		   so as to be between 100*ETA and 0.001.
		   Then the computed value h0 is given by

			  h0^(-2) = 1. / ( tol * w0^2 ) + tol * ( norm(f) )^2

		   where   w0     = max( fabs(*t), fabs(tout) ),
				   f      = the initial value of the vector f(t,y), and
				   norm() = the weighted vector norm used throughout, given by
							the vmnorm function routine, and weighted by the
							tolerances initially loaded into the ewt array.

		   The sign of h0 is inferred from the initial values of tout and *t.
		   fabs(h0) is made < fabs(tout-*t) in any case.
		*/
		if(h0 == 0.) {
			tdist = fabs(tout - *t);
			w0    = std::max(fabs(*t), fabs(tout));
			if(tdist < 2. * ETA * w0) {
				fprintf(stderr, "[lsoda] tout too close to t to start integration\n ");
				terminate(istate);
				return;
			}
			tol = rtol_[1];
			if(itol_ > 2) {
				for(size_t i = 2; i <= n; i++)
					tol = std::max(tol, rtol_[i]);
			}
			if(tol <= 0.) {
				atoli = atol_[1];
				for(size_t i = 1; i <= n; i++) {
					if(itol_ == 2 || itol_ == 4)
						atoli = atol_[i];
					ayi = fabs(y[i]);
					if(ayi != 0.)
						tol = std::max(tol, atoli / ayi);
				}
			}
			tol = std::max(tol, 100. * ETA);
			tol = std::min(tol, 0.001);
			sum = vmnorm(n, yh_[2], ewt);
			sum = 1. / (tol * w0 * w0) + tol * sum * sum;
			h0  = 1. / sqrt(sum);
			h0  = std::min(h0, tdist);
			h0  = h0 * ((tout - *t >= 0.) ? 1. : -1.);
		} /* end if ( h0 == 0. )   */
		/*
		   Adjust h0 if necessary to meet hmax bound.
		*/
		rh = fabs(h0) * hmxi;
		if(rh > 1.)
			h0 /= rh;

		/*
		   Load h_ with h0 and scale yh_[2] by h0.
		*/
		h_ = h0;
		for(size_t i = 1; i <= n; i++)
			yh_[2][i] *= h0;
	} /* if ( *istate == 1 )   */
	/*
	   Block d.
	   The next code block is for continuation calls only ( *istate = 2 or 3 )
	   and is to check stop conditions before taking a step.
	*/
	if(*istate == 2 || *istate == 3) {
		nslast = nst;
		switch(itask) {
			case 1:
				if((tn_ - tout) * h_ >= 0.) {
					intdy(tout, 0, y, &iflag);
					if(iflag != 0) {
						fprintf(stderr,
							"[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask,
							tout);
						terminate(istate);
						return;
					}
					*t      = tout;
					*istate = 2;
					illin   = 0;
					return;
				}
				break;
			case 2:
				break;
			case 3:
				tp = tn_ - hu * (1. + 100. * ETA);
				if((tp - tout) * h_ > 0.) {
					fprintf(
						stderr, "[lsoda] itask = %d and tout behind tcur - hu\n", itask);
					terminate(istate);
					return;
				}
				if((tn_ - tout) * h_ < 0.)
					break;
				successreturn(after_step, y, t, itask, ihit, tcrit, istate);
				return;
			case 4:
				tcrit = rworks[0];
				if((tn_ - tcrit) * h_ > 0.) {
					fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
					terminate(istate);
					return;
				}
				if((tcrit - tout) * h_ < 0.) {
					fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tout\n");
					terminate(istate);
					return;
				}
				if((tn_ - tout) * h_ >= 0.) {
					intdy(tout, 0, y, &iflag);
					if(iflag != 0) {
						fprintf(stderr,
							"[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask,
							tout);
						terminate(istate);
						return;
					}
					*t      = tout;
					*istate = 2;
					illin   = 0;
					after_step(*t, ++y.begin());
					return;
				}
				break;
			case 5:
				if(itask == 5) {
					tcrit = rworks[0];
					if((tn_ - tcrit) * h_ > 0.) {
						fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
						terminate(istate);
						return;
					}
				}
				hmx  = fabs(tn_) + fabs(h_);
				ihit = fabs(tn_ - tcrit) <= (100. * ETA * hmx);
				if(ihit) {
					*t = tcrit;
					successreturn(after_step, y, t, itask, ihit, tcrit, istate);
					return;
				}
				tnext = tn_ + h_ * (1. + 4. * ETA);
				if((tnext - tcrit) * h_ <= 0.)
					break;
				h_ = (tcrit - tn_) * (1. - 4. * ETA);
				if(*istate == 2)
					jstart = -2;
				break;
		} /* end switch   */
	}     /* end if ( *istate == 2 || *istate == 3 )   */
	/*
	   Block e.
	   The next block is normally executed for all calls and contains
	   the call to the one-step core integrator stoda.

	   This is a looping point for the integration steps.

	   First check for too many steps being taken, update ewt ( if not at
	   start of problem).  Check for too much accuracy being requested, and
	   check for h_ below the roundoff level in *t.
	*/
	while(1) {
		if(*istate != 1 || nst != 0) {
			if((nst - nslast) >= mxstep) {
				std::cerr << "[lsoda] " << mxstep << " steps taken before reaching tout"
					 << std::endl;
				*istate = -1;
				terminate2(y, t);
				return;
			}

			ewset(yh_[1]);
			for(size_t i = 1; i <= n; i++) {
				if(ewt[i] <= 0.) {
					std::cerr << "[lsoda] ewt[" << i << "] = " << ewt[i] << " <= 0." << std::endl;
					*istate = -6;
					terminate2(y, t);
					return;
				}
				ewt[i] = 1. / ewt[i];
			}
		}
		tolsf = ETA * vmnorm(n, yh_[1], ewt);
		if(tolsf > 0.01) {
			tolsf = tolsf * 200.;
			if(nst == 0) {
				fprintf(stderr, "lsoda -- at start of problem, too much accuracy\n");
				fprintf(stderr, "         requested for precision of machine,\n");
				fprintf(stderr, "         suggested scaling factor = %g\n", tolsf);
				terminate(istate);
				return;
			}
			fprintf(stderr, "lsoda -- at t = %g, too much accuracy requested\n", *t);
			fprintf(stderr, "         for precision of machine, suggested\n");
			fprintf(stderr, "         scaling factor = %g\n", tolsf);
			*istate = -2;
			terminate2(y, t);
			return;
		}

		if((tn_ + h_) == tn_) {
			nhnil++;
			if(nhnil <= mxhnil) {
				fprintf(stderr, "lsoda -- warning..internal t = %g and h_ = %g are\n",
					tn_, h_);
				fprintf(stderr,
					"         such that in the machine, t + h_ = t on the next step\n");
				fprintf(stderr, "         solver will continue anyway.\n");
				if(nhnil == mxhnil) {
					std::cerr << "lsoda -- above warning has been issued " << nhnil
						 << " times, " << std::endl
						 << "       it will not be issued again for this problem" << std::endl;
				}
			}
		}

		/* Call stoda */
		stoda(neq, y, derivs, _data);
		if(kflag == 0) {
			/*
			   Block f.
			   The following block handles the case of a successful return from the
			   core integrator ( kflag = 0 ).
			   If a method switch was just made, record tsw, reset maxord,
			   set jstart to -1 to signal stoda to complete the switch,
			   and do extra printing of data if ixpr = 1.
			   Then, in any case, check for stop conditions.
			*/
			if (tn_ < tout) after_step(tn_, ++y.begin());

			init = 1;
			if(meth_ != mused) {
				tsw    = tn_;
				maxord = mxordn;
				if(meth_ == 2)
					maxord = mxords;
				jstart = -1;
				if(ixpr) {
					if(meth_ == 2)
						std::cerr << "[lsoda] a switch to the stiff method has occurred "
							 << std::endl;
					if(meth_ == 1)
						std::cerr << "[lsoda] a switch to the nonstiff method has occurred"
							 << std::endl;
				}
			} /* end if ( meth_ != mused )   */
			/*
			   itask = 1.
			   If tout has been reached, interpolate.
			*/
			if(1 == itask) {
				if((tn_ - tout) * h_ < 0.)
					continue;

				intdy(tout, 0, y, &iflag);
				*t      = tout;
				*istate = 2;
				illin   = 0;
				after_step(*t, ++y.begin());
				return;
			}
			/*
			   itask = 2.
			*/
			if(itask == 2) {
				successreturn(after_step, y, t, itask, ihit, tcrit, istate);
				return;
			}
			/*
			   itask = 3.
			   Jump to exit if tout was reached.
			*/
			if(itask == 3) {
				if((tn_ - tout) * h_ >= 0.) {
					successreturn(after_step, y, t, itask, ihit, tcrit, istate);
					return;
				}
				continue;
			}
			/*
			   itask = 4.
			   See if tout or tcrit was reached.  Adjust h_ if necessary.
			*/
			if(itask == 4) {
				if((tn_ - tout) * h_ >= 0.) {
					intdy(tout, 0, y, &iflag);
					*t      = tout;
					*istate = 2;
					illin   = 0;
					after_step(*t, ++y.begin());
					return;
				}
				else {
					hmx  = fabs(tn_) + fabs(h_);
					ihit = fabs(tn_ - tcrit) <= (100. * ETA * hmx);
					if(ihit) {
						successreturn(after_step, y, t, itask, ihit, tcrit, istate);
						return;
					}
					tnext = tn_ + h_ * (1. + 4. * ETA);
					if((tnext - tcrit) * h_ <= 0.)
						continue;
					h_     = (tcrit - tn_) * (1. - 4. * ETA);
					jstart = -2;
					continue;
				}
			} /* end if ( itask == 4 )   */
			/*
			   itask = 5.
			   See if tcrit was reached and jump to exit.
			*/
			if(itask == 5) {
				hmx  = fabs(tn_) + fabs(h_);
				ihit = fabs(tn_ - tcrit) <= (100. * ETA * hmx);
				successreturn(after_step, y, t, itask, ihit, tcrit, istate);
				return;
			}

		} /* end if ( kflag == 0 )   */
		/*
		   kflag = -1, error test failed repeatedly or with fabs(h_) = hmin.
		   kflag = -2, convergence failed repeatedly or with fabs(h_) = hmin.
		*/
		if(kflag == -1 || kflag == -2) {
			fprintf(stderr, "lsoda -- at t = %g and step size h_ = %g, the\n", tn_, h_);
			if(kflag == -1) {
				fprintf(stderr, "         error test failed repeatedly or\n");
				fprintf(stderr, "         with fabs(h_) = hmin\n");
				*istate = -4;
			}
			if(kflag == -2) {
				fprintf(stderr, "         corrector convergence failed repeatedly or\n");
				fprintf(stderr, "         with fabs(h_) = hmin\n");
				*istate = -5;
			}
			big   = 0.;
			imxer = 1;
			for(size_t i = 1; i <= n; i++) {
				size = fabs(acor[i]) * ewt[i];
				if(big < size) {
					big   = size;
					imxer = i;
				}
			}
			terminate2(y, t);
			return;
		} /* end if ( kflag == -1 || kflag == -2 )   */
	}     /* end while   */
} /* end lsoda   */


/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Simpler interface.
 *
 * @Param f System
 * @Param neq, size of system.
 * @Param y, init values of size neq
 * @Param yout, results vector for size neq+1, ignore yout[0]
 * @Param t, start time.
 * @Param tout, stop time.
 * @Param _data
 * @Param rtol, relative tolerance.
 * @Param atol, absolute tolerance.
 */
/* ----------------------------------------------------------------------------*/
template <class Functor, class AfterStep>
void LSODA::lsoda_update(Functor &derivs, AfterStep &after_step, const size_t neq, std::vector<double> &y,
	/*std::vector<double> &yout_notused,*/ double *t, const double tout, /*int *istate,*/ void *_data,
	double rtol, double atol)
{
	std::array<int, 7> iworks    = {{0}};
	std::array<double, 4> rworks = {{0.0}};

	int itask, iopt, jt;

	// cout << "Debug : rtol " << rtol << ". atol " << atol << std::endl;

	itask = 1;
	iopt  = 0;
	jt    = 2;

	yout.resize(neq + 1);

	// Set the tolerance. We should do it only once.
	rtol_.resize(neq + 1, rtol);
	atol_.resize(neq + 1, atol);
	rtol_[0] = 0;
	atol_[0] = 0;

	// Fill-in values.
	for(size_t i = 1; i <= neq; i++)
		yout[i] = y[i - 1];

	lsoda(derivs, after_step, neq, yout, t, tout, itask, &iState, iopt, jt, iworks, rworks, _data);
	
	for (int i=0; i<y.size(); ++i) y[i] = yout[i+1];
}



