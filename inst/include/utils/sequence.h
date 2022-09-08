#ifndef UTILS_MATH_SEQUENCE_H_
#define UTILS_MATH_SEQUENCE_H_

#include <vector>
#include <cmath>

std::vector<double> my_log_seq(double x0, double xf, int N){
	std::vector<double> grid;
	for (int i=0; i<N; ++i) grid.push_back(exp(log(x0) + (double(i)/(N-1))*(log(xf)-log(x0))));
	return grid;
}

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

#endif

