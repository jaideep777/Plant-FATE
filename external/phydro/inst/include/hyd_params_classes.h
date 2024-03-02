#ifndef PHYDRO_PARAMS_CLASSES_H
#define PHYDRO_PARAMS_CLASSES_H

#include <iostream>

namespace phydro{

class ParPlant{
	public:
	double conductivity;
	double psi50;
	double b;

	double tchome = 25;

	double h_canopy = 20;
	double h_wind_measurement = 22;

	inline ParPlant(double _conductivity, double _psi50, double _b){
		conductivity = _conductivity;
		psi50 = _psi50;
		b = _b;
	}

	inline ParPlant(double _conductivity, double _psi50, double _b,
	         double _tchome, double _h_canopy, double _h_wind_measurement){

		conductivity = _conductivity;
		psi50 = _psi50;
		b = _b;

		tchome = _tchome;
		h_canopy = _h_canopy;
		h_wind_measurement = _h_wind_measurement;
	}

	inline void print(){
		std::cout << "ParPlant:\n";
		std::cout << "   conductivity = " << conductivity << '\n';
		std::cout << "   psi50 = " << psi50 << '\n';
		std::cout << "   b = " << b << '\n';
		std::cout << "   tchome = " << tchome << '\n';
		std::cout << "   h_canopy = " << h_canopy << '\n';
		std::cout << "   h_wind_measurement = " << h_wind_measurement << '\n';
	}

};


class ParCost{
	public:
	double alpha;
	double gamma;

	inline ParCost(double _a, double _g){
		alpha = _a;
		gamma = _g;
	}
};


} // phydro

#endif
