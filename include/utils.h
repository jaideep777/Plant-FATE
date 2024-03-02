#ifndef FLARE_UTILS_UTILS_H
#define FLARE_UTILS_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace flare{

// print a std::vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	//os << std::setprecision(12);
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	// os << '\n';
	return os;
}

namespace utils{

// https://stackoverflow.com/questions/14997165/fastest-way-to-get-a-positive-modulo-in-c-c
inline double positive_fmod(double x, double period){
	double res = fmod(x, period);
	if (res < 0) res += period;
	return res;
}

inline int positive_mod(int x, int period){
	int res = x % period;
	if (res < 0) res += period;
	return res;
}

inline std::string to_lower(std::string s){
	std::transform(s.begin(), s.end(), s.begin(),
			[](unsigned char c){ return std::tolower(c); });
	return s;
}


} // namespace utils
} // namespace flare

#endif
