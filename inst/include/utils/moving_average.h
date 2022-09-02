#ifndef UTILS_MATH_MOVING_AVG_H_
#define UTILS_MATH_MOVING_AVG_H_

#include <iostream>

#include <list>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

class MovingAverager{
	private:
	double area_sum = 0;
	std::list<double> areas;
	std::list<double> dts;
	std::list<double> f_hist;
	std::list<double> t_hist;

	public:
	// seed output history
	double T = 5;
	bool debug = true;

	public:
	
	inline void set_interval(double _T){
		T = _T;
	}
	
	inline void push(double t, double f){
		if (f_hist.empty()){
			f_hist.push_back(f);
			t_hist.push_back(t);
			//areas.push_back(0);
			return;
		}
		
		assert(t > t_hist.back());

		double f_lo = f_hist.back();
		double t_lo = t_hist.back();
		double area_new = (f + f_lo) / 2 * (t - t_lo);
		f_hist.push_back(f);
		t_hist.push_back(t);
		dts.push_back(t - t_lo);
		areas.push_back(area_new);
		area_sum += area_new;
		
		while(!t_hist.empty()){
			if (t_hist.front() < t - T){  // past t is beyond the averaging interval
				area_sum -= areas.front();
				f_hist.pop_front();
				t_hist.pop_front();
				areas.pop_front();
			}
			else{
				break;
			}
		}
	}

	inline double get(){
		if (t_hist.size() == 0) return 0;
		else if (t_hist.size() == 1) return f_hist.front();
		else return area_sum/(t_hist.back()-t_hist.front());
	}
	
	inline void clear(){
		area_sum = 0;
		f_hist.clear();
		t_hist.clear();
		areas.clear();
	}
	
	inline void print(){
		std::cout << "MovingAverager:\n";
		auto f_it = f_hist.begin();
		auto t_it = t_hist.begin();
		auto a_it = areas.begin();
		auto d_it = dts.begin();
		
		std::cout << "   t\tf\tA\tdt\n";
		while(f_it != f_hist.end()){
			std::cout << "   " << *t_it << "\t" << *f_it << "\t";
			std::cout << ((a_it != areas.end())? *a_it : 0) << "\t";
			std::cout << ((d_it != dts.end())? *d_it : 0) << "\n";
			++f_it; ++t_it; ++a_it;
		}
		std::cout << "-----\n";
		std::cout << "sum = " << area_sum << ", avg = " << get() << "\n";
		
	}

	inline void print_summary(){
		std::cout << "MovingAverager:  t = " << t_hist.front() << " " << t_hist.back() << ", x = " << f_hist.front() << " " << f_hist.back() << ", xmean = " << get() << ", (" << t_hist.back() - t_hist.front() << ")\n";
	}

	inline double get_exp(double c = 0){
		if (t_hist.size() == 0) return 0;
		if (t_hist.size() == 1) return f_hist.front();

		if (debug) std::cout << "Exp avg: ";
		auto f_it = f_hist.rbegin();
		auto t_it = t_hist.rbegin();
		double avg = 0;
		double D = 0;
		double t0 = *t_hist.rbegin();
		int i=0;
		while (f_it != prev(f_hist.rend())){
			double t_hi = *t_it;
			double t_lo = *next(t_it);
			double f_hi = (*f_it);
			double w_hi = exp(-c*(t0-t_hi));
			double f_lo = (*next(f_it));
			double w_lo = exp(-c*(t0-t_lo));
			double dt = t_hi - t_lo;
			avg += (f_hi*w_hi + f_lo*w_lo)/2 * dt;
			D += (w_hi + w_lo)/2 * dt;
			++f_it; ++t_it;

			if (debug) std::cout << (f_hi*w_hi + f_lo*w_lo)/2 * dt << " ";
		}
		avg /= D;

		if (debug) std::cout << " | " << avg << std::endl;

		return avg;
	}


	inline double get_cesaro(int nmax, double c=0){
		if (t_hist.size() == 0) return 0;
		if (t_hist.size() == 1) return f_hist.front();

		// f_hist has 2 or more elements, so compute
		nmax = std::min(int(nmax), int(f_hist.size()));
		std::vector<double> avgs(nmax, 0);

		// assert that t is equally spaced
		auto tt = std::minmax_element(dts.begin(), dts.end());
		double thresh = 0.01;
		if (abs(*tt.first - *tt.second) > thresh) throw std::runtime_error("Cesaro Average requires equally spaced intervals.");

		for (int n=1; n<nmax; ++n){
			avgs[n] = 0;
			double D = 0;
			auto f_it = f_hist.rbegin();
			for (int i=0; i<n; ++i, ++f_it){
				avgs[n] += (*f_it) * exp(-c*i);
				D += exp(-c*i);
			}
			avgs[n] /= D;
		}

		double avg = 0;
		double D = 0;
		for (int i=1; i<nmax; ++i){
			avg += avgs[i]*exp(-c*(i-1)); // Note avgs[0] is not defined, the first number is at index 1
			D += exp(-c*i);
		}
		avg /= D;

		if (debug){
			std::cout << "Cesaro: ";
			for (auto x : avgs) std::cout << x << " ";
			std::cout << " | " << avg << std::endl;
		}

		return avg;
	}

};


#endif
