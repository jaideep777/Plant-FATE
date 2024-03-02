// Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
// adjust stepsize. Input are
//  y[...], dydx[...]  -- the dependent variable vector and its derivative at the starting value of the independent variable x.
//  htry, hdid         -- the stepsize to be attempted htry, and the stepsize that was actually accomplished hdid
//  hnext              -- the estimated next stepsize,
//  derivs             -- the user-supplied routine that computes the right-hand side derivatives.
//
// Also used in the class are the required accuracy eps, and the vector yscal[...] against which the error is scaled.
// On output, y and x are replaced by their new values
template <class functor>
void RKCK45::RKStep(container& y, container& dydx, double& x, double htry,
                    double& hdid, double& hnext, functor& derivs){

	container yerr(y.size()), ytemp(y.size());	// TODO: Should these also be made class members to prevent reallocation?
	double errmax;
	double h=htry; // Set stepsize to the initial trial value.
	for (;;) {     // infinite loop
		RKTry(y, dydx, x, h, ytemp, yerr, derivs);// Take a step.

		errmax=0.0;                               // Evaluate accuracy.
		for (int i=0; i<y.size(); i++) errmax=std::max(errmax, fabs(yerr[i]/yscal[i]));
		//errmax /= eps;                            // Scale relative to required tolerance.
		// std:: cout << "x = " << x << ", h = " << h << ", errmax = " << errmax << "\n";
		if (errmax <= 1.1) break;                 // Step succeeded. Compute size of next step.
		double htemp=h*fmax(SAFETY*pow(errmax,PSHRNK), 0.2); // Truncation error too large, reduce stepsize, max by a factor of 0.2.
		h=htemp; //(h >= 0.0 ? std::max(htemp,0.1*h) : std::min(htemp,0.1*h)); // No more than a factor of 10.
		if ((x+h) == x) std::cerr<<"stepsize underflow in rkqs"<<std::endl;
	}
	if (errmax < 0.5) hnext=h*fmin(SAFETY*pow(errmax,PGROW), 5); // Step is too small, increase it next time, no more than 5 times
	else hnext=h;                                      
	// std::cout << "hnext = " << hnext << "\n";
	hdid=h;
	x += h;
	for (int i=0; i<y.size();i++) y[i] = ytemp[i];
}

// Given values for variables y[...] and their derivatives dydx[...] known at x, use
// the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h
// and return the incremented variables as yout[..]. Also return an estimate of the local
// truncation error in yout using the embedded fourth-order method. The user supplies the routine
// derivs(x,y,dydx), which returns derivatives dydx at x.
template <class functor>
void RKCK45::RKTry(container& y, container& dydx, double& x, double h, container& yout, container& yerr, functor& derivs){ 
	//std::cout << "try:\n";
	//                              a0 a1   a2   a3   a4   a5
	static constexpr double as[] = {0, 0.2, 0.3, 0.6, 1.0, 0.875};
	//                              c0          c1   c2           c3           c4   c5
	static constexpr double cs[] = {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0};
	//                              dc0                   dc1  dc2                    dc3                     dc4              dc5
	static constexpr double dc[] = {cs[0]-2825.0/27648.0, 0.0, cs[2]-18575.0/48384.0, cs[3]-13525.0/55296.0, -277.00/14336.0, cs[5]-0.25};
	static constexpr double bs[6][6] =
		{{0.0,            0.0,         0.0,           0.0,              0.0,          0.0},
		 {0.2,            0.0,         0.0,           0.0,              0.0,          0.0},
		 {3.0/40.0,       9.0/40.0,    0.0,           0.0,              0.0,          0.0},
		 {0.3,           -0.9,         1.2,           0.0,              0.0,          0.0},
		 {-11.0/54.0,     2.5,        -70.0/27.0,     35.0/27.0,        0.0,          0.0},
		 {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, 0.0}};
  
	int N=y.size();
	
	container& k0 = dydx; // just different name for the same variable
	
	// First step.	
	for (int i=0; i<N; i++) yt[i] = y[i] + h*bs[1][0]*k0[i];
	derivs(x+as[1]*h, yt.begin(), k1.begin(), nullptr);                            
	++nfe;

	// Second step.
	for (int i=0; i<N; i++) yt[i] = y[i] + h*(bs[2][0]*k0[i]+bs[2][1]*k1[i]);
	derivs(x+as[2]*h,yt.begin(),k2.begin(), nullptr);                             
	++nfe;

	// Third step.
	for (int i=0; i<N; i++) yt[i] = y[i] + h*(bs[3][0]*k0[i]+bs[3][1]*k1[i]+bs[3][2]*k2[i]);
	derivs(x+as[3]*h,yt.begin(),k3.begin(), nullptr);                             
	++nfe;	

	// Fourth step.
	for (int i=0; i<N; i++) yt[i] = y[i] + h*(bs[4][0]*k0[i]+bs[4][1]*k1[i]+bs[4][2]*k2[i]+bs[4][3]*k3[i]);
	derivs(x+as[4]*h,yt.begin(),k4.begin(), nullptr);                             
	++nfe;	

	// Fifth step.
	for (int i=0;i<N;i++)   yt[i] = y[i] + h*(bs[5][0]*k0[i]+bs[5][1]*k1[i]+bs[5][2]*k2[i]+bs[5][3]*k3[i]+bs[5][4]*k4[i]);
	derivs(x+as[5]*h,yt.begin(),k5.begin(), nullptr);                             
	++nfe;	

	// Sixth step.
	// Accumulate increments with proper weights.
	for (int i=0;i<N;i++)   yout[i] = y[i] + h*(cs[0]*k0[i]+cs[2]*k2[i]+cs[3]*k3[i]+cs[5]*k5[i]);
	
	// Estimate error as difference between fourth and fifth order methods.
	for (int i=0; i<N; i++) yerr[i] = h*(dc[0]*dydx[i]+dc[2]*k2[i]+dc[3]*k3[i]+dc[4]*k4[i]+dc[5]*k5[i]);
}

// The main function which takes next RK step with predefined accuracy
// returns bool:   true if time becomes larger of equal to t_stop, otherwise false
// arguments:
//           1) x      -- current time (independent variable), which will be changes for a step when successful
//           2) y      -- current solution (dependent variable)
//           3) derivs -- function which evaluates derivatives
template <class functor>
void RKCK45::Step(double& x, container& y, functor& derivs, double hmax){
	
	// resize the container if necessary
	if (y.size() != sys_size) resize(y.size());

	// Calculates derivatives for the new step
	//cout << "\n>> init derivs ~~~\n";
	derivs(xt,y.begin(),dydx.begin(), nullptr);
	++nfe;
	// good way of determining desired accuracy
	for (int i=0; i<y.size(); i++) yscal[i] =  eps_abs + eps_rel*(a_y * fabs(y[i]) + a_dydt * fabs(dydx[i]*ht)); // fabs(y[i]) + fabs(dydx[i]*ht) + 1e-3;
	// std::cout << "eps = " << eps_abs << "/" << eps_rel << "\n";
	double hnext, hdid;
	// Cals the Runge-Kutta rountine
	// takes: y -- dependent variable, dydx -- derivative at the beginning
	//        ht -- trial setpsize, hdid -- 
	ht = std::min(ht, hmax);
	ht = std::max(ht, hmin);
	RKStep(y, dydx, xt, ht, hdid, hnext, derivs); // Make one RK step
	if (hdid == ht) ++nok; else ++nbad; // good or bad step
	ht = hnext;
	x = xt;
	//return x<t_stop;
}
 
	//template <class functor>
	//void Step_RK4(double &x, double h, container& y, functor& derivs){
		////static container k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());// temporary arrays
		////static container yt(y.size());
		//// resize the container if necessary
		//if (y.size() != sys_size) resize(y.size());

		//double h2=h*0.5;
		//double xh = x + h2;
		//derivs(x, y, k1);     // First step : evaluating k1
		//for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k1[i];// Preparing second step by  ty <- y + k1/2
		//derivs(xh, yt, k2);                                    // Second step : evaluating k2
		//for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k2[i];// Preparing third step by   yt <- y + k2/2
		//derivs(xh, yt, k3);                                    // Third step : evaluating k3
		//for (int i=0; i<y.size(); i++) yt[i] = y[i] +  h*k3[i];// Preparing fourth step  yt <- y + k3
		//derivs(x+h, yt, k4);                                   // Final step : evaluating k4
		//for (int i=0; i<y.size(); i++) y[i] += h/6.0*(k1[i]+2.0*(k2[i]+k3[i])+k4[i]);

		//x += h;
	//}

template <class functor, class AfterStep>
void RKCK45::Step_to(double t_stop, double& x, container& y, functor& derivs, AfterStep &after_step){
	// std::cout << "ode step: " << x << " --> " << t_stop << '\n';
	while (x < t_stop){
		// Advance x without stepping if t_stop is within hmin
		// this case mostly arises due to float comparisons in manual stepping
		// doing this prevents reset of h to hmin
		if ((t_stop - x) < hmin){
			x = t_stop; 
			break;
		}

		// Otherwise, take RK steps until tstop
		Step(x,y,derivs, t_stop-x);
		after_step(x, y.begin());
	}
}





