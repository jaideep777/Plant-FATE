#include <iostream>
#include <cmath>
#include "mcmc.h"
using namespace std;

// compile: g++ -o 1 mcmc_test.cpp

// Test with R:
// ------------------
// library(tidyverse)
// dat = read.delim("/home/jjoshi/codes/libpspm/tests/mcmc_test.txt", header=F)
// dat %>% ggplot(aes(x=V2, y=V3, col=V1)) + 
//   geom_point(alpha=0.5) + 
//   # geom_line(alpha=0.1) +
//   geom_point(x=1,y=2,col="red")+
//   geom_point(x=5,y=2,col="red")+
//   geom_point(x=1,y=10,col="red")+
//   geom_point(x=5,y=10,col="red")

// hist(dat$V2[-(1:200)], breaks=50)
// hist(dat$V3[-(1:200)], breaks=50)


int main() {

	auto targetDistribution = [](const std::vector<double>& x) {
        std::vector<double> means = {1,2};
		double result = 1.0;
        for (int i=0; i<x.size(); ++i) {
            result *= std::exp(-(x[i]-means[i])* (x[i]-means[i])) + std::exp(-(x[i]-5*means[i])* (x[i]-5*means[i])); // Example: Multivariate Normal Distribution
        }
        return result;
    };

	auto targetDistribution_gaussian = [](const std::vector<double>& x) {
        std::vector<double> means = {1,2};
		std::vector<std::vector<double>> sd = {{1, 0.5}, {0.5, 1}};
		double cov_det = sd[0][0] * sd[1][1] - sd[1][0] * sd[0][1];
		std::vector<std::vector<double>> sd_inv = {{sd[1][1]/cov_det, -sd[1][0]/cov_det}, {-sd[0][1]/cov_det , sd[0][0]/cov_det}};
		

		double a = x[0] - means[0];
		double b = x[1] - means[1];

		double mah_distance = a * a *sd_inv[0][0] + a * b *(sd_inv[0][1] + sd_inv[1][0]) + b * b *sd_inv[1][1];

		double result = std::exp(-0.5 * mah_distance) / sqrt(pow(2*M_PI, 2) * cov_det);
        return result;
    };

	std::cout << "Creating sampler" << std::endl;

	MCMCSampler sampler({-5,-5}, {5,5}, {0.25, 0.25}, 4, 1000, 1);

	std::cout << "Created sampler" << std::endl;
	std::cout << "Running chains" << std::endl;
	sampler.run_chains(targetDistribution_gaussian, 10000);
	std::cout << "Chains ran" << std::endl;

	std::cout << "Writing output" << std::endl;

	std::ofstream fout("mcmc_test.txt");
	sampler.chainList[0].printSamples(fout);
	sampler.chainList[1].printSamples(fout);
	sampler.chainList[2].printSamples(fout);
	sampler.chainList[3].printSamples(fout);
	fout.close();

	std::cout << "Output created" << std::endl;
	std::cout << "Running gelman-rubin test" << std::endl;

	std::vector<double> rubintest = sampler.gelman_rubin_test();

	std::cout << "GELMAN RUBIN TEST" << std::endl;
	std::cout <<  rubintest[0] << "\t"  << rubintest[1] << std::endl;

	return 0;
}
