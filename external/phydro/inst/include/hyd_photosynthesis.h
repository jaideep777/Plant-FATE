#ifndef PHYDRO_PHOTOSYNTHESIS_H
#define PHYDRO_PHOTOSYNTHESIS_H

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "temperature_dependencies_photosynthesis.h"

namespace phydro{



class ParPhotosynth{
	public:
	double kmm;
	double gammastar;
	double phi0;
	double ca;     // Partial pressure of CO2 [Pa]
	double delta;  // TODO: Replace name with brd / rdark

	FtempVcmaxJmaxMethod ftemp_vj_method;
	FtempRdMethod        ftemp_rd_method;
	FtempBrMethod        ftemp_br_method;

	double Iabs;  // Net absorbed PAR [umol m-2 s-1]
	double patm;  // Atmospheric pressure [Pa]

	double fT_vcmax;
	double fT_jmax;
	double fT_rd;

	inline ParPhotosynth(double _tc, double _patm, double _kphio, double _co2, double _ppfd, double _fapar, double _rdark25, double _tcgrowth, double _tchome,
				  FtempVcmaxJmaxMethod _ftemp_vj_method = FV_kumarathunge19, 
				  FtempRdMethod        _ftemp_rd_method = FR_heskel16, 
				  FtempBrMethod        _ftemp_br_method = FB_atkin15){
		
		ftemp_vj_method = _ftemp_vj_method;
		ftemp_rd_method = _ftemp_rd_method;
		ftemp_br_method = _ftemp_br_method;

		fT_vcmax = calc_ftemp_inst_vcmax(_tc, _tcgrowth, 25.0, ftemp_vj_method);
		fT_jmax  = calc_ftemp_inst_jmax(_tc, _tcgrowth, _tchome, 25, ftemp_vj_method);
		fT_rd    = calc_ftemp_inst_rd(_tc, _ftemp_rd_method);

		kmm = calc_kmm(_tc, _patm);
		gammastar = calc_gammastar(_tc, _patm);

		phi0 = _kphio*calc_ftemp_kphio(_tc);
		Iabs = _ppfd * _fapar;
		ca = _co2 * _patm * 1e-6;
		patm = _patm;
		delta = _rdark25 * fT_rd / fT_vcmax;

	}

	inline void print(){
		std::cout << "Par Photosynth:\n";
		std::cout << "   fT_vcmax = " << fT_vcmax << '\n';
		std::cout << "   fT_jmax = " << fT_jmax << '\n';
		std::cout << "   fT_rd = " << fT_rd << '\n';
		std::cout << "   kmm = " << kmm << '\n';
		std::cout << "   gammastar = " << gammastar << '\n';
		std::cout << "   phi0 = " << phi0 << '\n';
		std::cout << "   Iabs = " << Iabs << '\n';
		std::cout << "   ca = " << ca << '\n';
		std::cout << "   patm = " << patm << '\n';
		std::cout << "   delta = " << delta << '\n';
		std::cout << "   ftemp_vj_method = " << ftemp_vj_method << '\n';
		std::cout << "   ftemp_rd_method = " << ftemp_rd_method << '\n';
		std::cout << "   ftemp_br_method = " << ftemp_br_method << '\n';
	}
};



struct ACi{
	double a;
	double ci;
	bool isVcmaxLimited;
};


inline double QUADM(double A, double B, double C){
	if (A==0) return -C/B;
	else return (-B - sqrt(B*B - 4*A*C)) / (2*A);
}

inline double QUADP(double A, double B, double C){
	if (A==0) return -C/B;
	else return (-B + sqrt(B*B - 4*A*C)) / (2*A);
}


inline ACi calc_assim_rubisco_limited(double _gs, double vcmax, ParPhotosynth par_photosynth){
	double ca = par_photosynth.ca;            // ca is in Pa
	double gs = _gs * 1e6/par_photosynth.patm; // convert to umol/m2/s/Pa

	double d = par_photosynth.delta;

	double A = -1.0 * gs;
	double B = gs * ca - gs * par_photosynth.kmm - vcmax*(1-d);
	double C = gs * ca * par_photosynth.kmm + vcmax * (par_photosynth.gammastar + par_photosynth.kmm*d);

	ACi res;
	res.ci = QUADM(A,B,C);
	res.a  = gs*(ca-res.ci);
	res.isVcmaxLimited = true;

	return res;

}

// FIXME: rewrite in terms of J and Jmax
inline ACi calc_assim_light_limited(double _gs, double jmax, ParPhotosynth par_photosynth){
	double ca = par_photosynth.ca;             // ca is in Pa
	double gs = _gs * 1e6/par_photosynth.patm;  // convert to umol/m2/s/Pa

	//gs += 1e-12; // FIXME: delete?

	double d = par_photosynth.delta;

	double phi0iabs = par_photosynth.phi0 * par_photosynth.Iabs;
	double jj = 4*phi0iabs/jmax;
	double jlim = phi0iabs / sqrt(1+ jj*jj);

	double A = -1.0 * gs;
	double B = gs * ca - gs * 2 * par_photosynth.gammastar - jlim*(1-d);
	double C = gs * ca * 2*par_photosynth.gammastar + jlim * (par_photosynth.gammastar + d*par_photosynth.kmm);

	ACi res;
	res.ci = QUADM(A,B,C);
	res.a  = gs*(ca-res.ci);
	res.isVcmaxLimited = false;

	return res;

}


// I have tried using a smooth minimum for calculating the limiting rate, but it absolutely does not work
// with the optimizer. I have no idea why.
inline ACi calc_assimilation_limiting(double vcmax, double jmax, double gs, ParPhotosynth par_photosynth){
	auto Ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth);
	auto Aj = calc_assim_light_limited(gs, jmax, par_photosynth);

	if (Ac.ci > Aj.ci ) return Ac; 
	else				return Aj;
}

} // phydro

#endif

