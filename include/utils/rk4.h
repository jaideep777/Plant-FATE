template <class functor, class container>
void Euler(double x, double h, container& y, functor& derivs){
	container fk(y.size());
	derivs(x, y, fk);
	for (int i=0; i<y.size(); i++) y[i] += h*fk[i]; 
}

template <class functor, class container>
void RK4(double x, double h, container& y, functor& derivs){
	static container k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());// temporary arrays
	static container yt(y.size());

	double h2=h*0.5;
	double xh = x + h2;
	derivs(x, y, k1);     // First step : evaluating k1
	for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k1[i];// Preparing second step by  ty <- y + k1/2
	derivs(xh, yt, k2);                                    // Second step : evaluating k2
	for (int i=0; i<y.size(); i++) yt[i] = y[i] + h2*k2[i];// Preparing third step by   yt <- y + k2/2
	derivs(xh, yt, k3);                                    // Third step : evaluating k3
	for (int i=0; i<y.size(); i++) yt[i] = y[i] +  h*k3[i];// Preparing fourth step  yt <- y + k3
	derivs(x+h, yt, k4);                                   // Final step : evaluating k4
	for (int i=0; i<y.size(); i++) y[i] += h/6.0*(k1[i]+2.0*(k2[i]+k3[i])+k4[i]);
}


