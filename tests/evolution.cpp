#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

class Species{
	public:
	bool isResident;

	double x,y; // traits
	double u;
	double t;
	double r = 1;
	double K;

	double invasion_fitness;
	double dlnudt;
	vector<Species*> probes;
	vector<double> fitness_gradient;
	double nudge_dx = 0.01;

	void set_traits(vector<double> tvec){
		x = tvec[0];
		y = tvec[1];
		K = fmax(1e-20, 1 - (x*x + y*y));
	}
	vector<double> get_traits(){
		return {x,y};
	}
	Species(vector<double> tvec, double u0, bool res){
		set_traits(tvec);
		fitness_gradient.resize(tvec.size(), 0);
		u = u0;
		isResident = res;
	}
	
	void nudgeTraits(){
		vector<double> traits_res = get_traits();
		vector<double> dx(traits_res.size());
		for (int i=0; i<dx.size(); ++i) dx[i] = fitness_gradient[i]*nudge_dx;

		for (int i=0; i<dx.size(); ++i) traits_res[i] += dx[i];
		set_traits(traits_res);

		for (auto m : probes){
			vector<double> traits_mut = m->get_traits();
			for (int j=0; j<traits_mut.size(); ++j) traits_mut[j] += dx[j];
			m->set_traits(traits_mut);
		}
	}

};

class Solver{
	public:
	double _dt = 0.01;
	double c = 0.18;
	double fg_dx = 0.1;

	vector<Species*> species_vec;
	double t =0;

	void createRandomResidents(int n){
		for (int i=0; i<n; ++i){
			double r = double(rand())/RAND_MAX;
			double theta = double(rand())/RAND_MAX*2*3.14159265;
			double x = r*cos(theta);
			double y = r*sin(theta);
			Species* s = new Species({x,y}, 0.02, true);
			species_vec.push_back(s);
		}
	}

	void createProbes(){
		for (auto spp : species_vec){
			vector<double> traits = spp->get_traits();
			for (int i=0; i<traits.size(); ++i){
				vector<double> traits_mutant = traits;
				traits_mutant[i] += fg_dx;
				Species * s = new Species(traits_mutant, 0.02, false);
				spp->probes.push_back(s);
			}
		}
		vector<Species*> species_vec_copy = species_vec;
		for (auto spp : species_vec_copy){
			for (auto m : spp->probes) species_vec.push_back(m);
		}
	}

	void reset(){
		t = 0;
	}

	// logistic model
	void step_to(double tf){
		while (t<tf){
			double dt = min(_dt, tf-t);

			// calc growth rates for each species 
			for (auto spp : species_vec){
				// cout << "spp = " << spp << "\n";
				double C = 0;
				for (auto s2 : species_vec){
					// no competition from mutants (including self)
					if (s2->isResident){
						double dx = (spp->x - s2->x);
						double dy = (spp->y - s2->y);
						double d = sqrt(dx*dx + dy*dy)/c;
						double a = exp(-d*d*d);

						C += a*s2->u;
						//cout << "   s2 = " << s2 << ", u = " << s2->u << ", K2 = " << s2->K << ", a = " << a << ", C_cumm = " << C << "\n";
					}
				}

				spp->dlnudt = spp->r*(1-C/spp->K);
				spp->invasion_fitness = spp->dlnudt;
			}

			// update density of each species
			for (auto spp : species_vec) spp->u = spp->u*exp(spp->dlnudt*dt);
			
			// update time
			t += dt;
		}

	}
};

int main(){
	
	Solver S;
	S.createRandomResidents(15);
	S.createProbes();

	// // add mutant probes
	// for (int i=0; i<1000; ++i){
	// 	double r = double(rand())/RAND_MAX;
	// 	double theta = double(rand())/RAND_MAX*2*3.14159265;
	// 	double x = r*cos(theta);
	// 	double y = r*sin(theta);
	// 	Species* s = new Species({x,y}, 0.02, false);
	// 	S.species_vec.push_back(s);
	// }

	ofstream ftrait("evol_traits.txt");
	ftrait << "iter\tID\tx\ty\tK\tu\tresident\tinvasion_fitness\tdfx\tdfy\t\n";

	ofstream fout("evolution.txt");
	fout << "iter\tt\t";
	for (int i=0; i<S.species_vec.size(); ++i) fout << "u" << i << "\t";
	fout << "\n";

	for (int iter=0; iter<100; ++iter){
		S.reset();

		for (double t=0; t<=20; t+=0.1){
			S.step_to(t);

			fout << iter << "\t" << t << "\t";
			for (int i=0; i<S.species_vec.size(); ++i) if (S.species_vec[i]->isResident) fout << S.species_vec[i]->u << "\t";
			fout << "\n";

		}

		for (auto spp : S.species_vec){
			spp->fitness_gradient.clear();
			for (auto m : spp->probes){
				spp->fitness_gradient.push_back((m->invasion_fitness - spp->invasion_fitness)/S.fg_dx);
			}
		}

		for (auto spp : S.species_vec) spp->nudgeTraits();

		for (int i=0; i<S.species_vec.size(); ++i){
			auto spp = S.species_vec[i];
			ftrait << iter << "\t" 
				<< i << "\t"
				<< spp->x << "\t"
				<< spp->y << "\t"
				<< spp->K << "\t"
				<< spp->u << "\t"
				<< spp->isResident << "\t"
				<< spp->invasion_fitness << "\t";
			if (spp->isResident){
				ftrait << spp->fitness_gradient[0] << "\t"
					<< spp->fitness_gradient[1] << "\t";
			}
			else{
				ftrait << "NA\tNA\t";
			}
			ftrait << "\n";

		}

	} 

	fout.close();
	ftrait.close();


}

