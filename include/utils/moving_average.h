#ifndef UTILS_MATH_MOVING_AVG_H_
#define UTILS_MATH_MOVING_AVG_H_

#include <iostream>

#include <list>
#include <cmath>

class MovingAverager{
	private:
	double area_sum = 0;
	std::list<double> areas;
	std::list<double> f_hist;
	std::list<double> t_hist;

	public:
	// seed output history
	double T = 5;
	
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
		
		while(f_it != f_hist.end()){
			std::cout << "   " << *t_it << "\t" << *f_it << "\t";
			std::cout << ((a_it != areas.end())? *a_it : 0) << "\n";
			++f_it; ++t_it; ++a_it;
		}
		std::cout << "-----\n";
		std::cout << "sum = " << area_sum << ", avg = " << get() << "\n";
		
	}

	inline void print_summary(){
		std::cout << "MovingAverager:  t = " << t_hist.front() << " " << t_hist.back() << ", x = " << get() << ", (" << t_hist.back() - t_hist.front() << ")\n";
	}

};


#endif
