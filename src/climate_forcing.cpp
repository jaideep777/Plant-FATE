#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "util.h"
#include "climate_forcing.h"

using namespace std;

namespace env{


int ClimateForcing::init(){
	std::string line;

	if (update_met){
		std::ifstream fin_met(metFile.c_str());
		if (!fin_met){
			throw std::runtime_error("Could not open file " + metFile);
		}

		// skip header
		getline(fin_met, line);

		while (fin_met.peek() != EOF){
			Clim clim1(clim);  // use copy constructor so that any user-defined values (e.g., co2) are preserved
			double t1;
			readNextLine_met(clim1, t1, fin_met);
			t_met.push_back(t1);
			v_met.push_back(clim1);
		}
		
		int n = t_met.size();
		if (n == 0) throw std::runtime_error("Climate time vector is empty");
		delta = t_met[n-1] - t_met[0] + (t_met[n-1]-t_met[0])/(n-1)*(1+1e-20); //+ 1e-20;
		deltaT = findMinDeltaT(t_met);
	}

	if (update_co2){
		std::ifstream fin_co2(co2File.c_str());
		if (!fin_co2){
			throw std::runtime_error("Could not open file " + co2File);
		}
		
		// skip header
		getline(fin_co2, line);
		
		// read CO2 file
		while (fin_co2.peek() != EOF){
			double co21;  // use copy constructor so that any user-defined values (e.g., co2) are preserved
			double t1;
			readNextLine_co2(co21, t1, fin_co2);
			
			t_co2.push_back(t1);
			v_co2.push_back(co21);
		}


	}
	
	return 0;

}


void ClimateForcing::set(double _tc, double _ppfd_max, double _ppfd, double _vpd, double _co2, double _elv, double _swp){
	clim.tc = _tc; 
	clim.ppfd_max = _ppfd_max; 
	clim.ppfd = _ppfd;
	clim.vpd = _vpd;
	clim.co2 = _co2;
	clim.elv = _elv;
	clim.swp = _swp;
};

void ClimateForcing::update_tc(double t, double _tc){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].tc = _tc;     
    }	
}

void ClimateForcing::update_ppfd_max(double t, double _ppfd_max){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].ppfd_max = _ppfd_max;     
    }	
}

void ClimateForcing::update_ppfd(double t, double _ppfd){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].ppfd = _ppfd;     
    }	
}

void ClimateForcing::update_vpd(double t, double _vpd){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].vpd = _vpd;     
    }
}

// void ClimateForcing::update_co2(double t, double _co2){
// 	int idx = idx_of(t);
//     if(idx >= 0){
//         v_met[idx].co2 = _co2;     
//     }
// }

void ClimateForcing::update_elv(double t, double _elv){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].elv = _elv;     
    }
	
}

void ClimateForcing::update_swp(double t, double _swp){
	int idx = idx_of(t);
    if(idx >= 0){
        v_met[idx].swp = _swp;    
    }
}


Clim ClimateForcing::interp(Clim &clim_prev, Clim &clim_next){
	return clim_prev;
}

int ClimateForcing::id(double t){
	int id = (t - *t_met.begin())/delta*t_met.size();
	id = id % t_met.size();
	return id;
}

int ClimateForcing::idx_of(double t){
	double tadj = std::fmod((t - *t_met.begin()),(*t_met.end() - *t_met.begin()));
    int i = floor((tadj - *t_met.begin())/deltaT);
    while((i+1 != t_met.size()) && t_met[i+1] < tadj) { ++i; }
    if(i == t_met.size()){
        return -1;
    }
    else{
        return i;
    }
}


void ClimateForcing::updateClimate(double t){

	// double tadj = t;  // adjusted t to lie between limits of observed data
	// while(tadj < *t_met.begin()) tadj += delta;
		
	if (update_met){
		int n = t_met.size();
		if (n == 1){
			clim = v_met[0];
		}
		else{
			double tadj = t;  // adjusted t to lie between limits of observed data
			while(tadj < t_met[0])   tadj += deltaT;
			int idx_now = idx_of(tadj);
			int idx_next = (idx_now+1) % t_met.size();
			clim = interp(v_met[idx_now], v_met[idx_next]);
		}
	}

	if (update_co2){
		int year = int(t);
		//std::cout << "CO2: " << t << " " << *t_co2.begin() << " " << *t_co2.rbegin() << "\n";
		if (year >=  *t_co2.begin() && year <= *t_co2.rbegin()){
			clim.co2 = v_co2[year - t_co2[0]];
			//std::cout << "Setting co2 @ t = " << year << " from (" << year - t_co2[0] << " / " << v_co2[year - t_co2[0]] << ")\n";
		}
		else if (year > *t_co2.rbegin()) {
			//std::cout << "Setting co2 @ end \n";
			clim.co2 = *v_co2.rbegin();
		}
	}
}



int ClimateForcing::readNextLine_met(Clim &clim, double &t, std::ifstream& fin_met){

	std::string                line, cell;
	
	// READ line
	std::getline(fin_met, line);

	std::stringstream          lineStream(line);
	
	// PARSE line
	std::getline(lineStream, cell, ',');
	int time = as<int>(cell);

	std::getline(lineStream, cell, ',');
	clim.tc = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim.vpd = as<double>(cell)*100;  // convert hPa to Pa
	
	std::getline(lineStream, cell, ',');
	clim.ppfd = as<double>(cell);

	std::getline(lineStream, cell, ',');
	clim.ppfd_max = as<double>(cell);
	
	std::getline(lineStream, cell, ',');
	clim.swp = -as<double>(cell);  // convert to negative (in file it is absolute value)

	t = time;
		
	return 0;
}

int ClimateForcing::readNextLine_co2(double &co2, double &t, std::ifstream& fin_co2){

	std::string                line, cell;
	
	// READ line
	std::getline(fin_co2, line);

	std::stringstream          lineStream(line);
	

	std::getline(lineStream, cell, ',');
	int time = as<double>(cell);

	std::getline(lineStream, cell, ',');
	double co2_ = as<double>(cell);
			

	// PARSE line
	t = time;
    co2 = co2_;
		
	return 0;
}

void ClimateForcing::print(){
	cout << "Climate meta:\n";
	cout << "  met_file = " << metFile << '\n';
	cout << "  co2_file = " << co2File << '\n';
	cout << "  met data will be updated from file: " << ((update_met)? "yes" : "no") << '\n';
	cout << "  CO2 will be updated from file: " << ((update_co2)? "yes" : "no") << '\n';
	cout << "  Current climate:\n";
	cout << "     tc      = " << clim.tc  << '\n';
	cout << "     ppfdmax = " << clim.ppfd_max << '\n';
	cout << "     ppfd    = " << clim.ppfd << '\n';
	cout << "     vpd     = " << clim.vpd << '\n';
	cout << "     co2     = " << clim.co2 << '\n';
	cout << "     elv     = " << clim.elv << '\n';
	cout << "     swp     = " << clim.swp << '\n';
}

void ClimateForcing::print_line(double t){
	int year = int(t);
	double month = (t-int(t))*12;
	std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
	std::cout << " | T/Im/I/D/CO2/Z/Ps = " << clim.tc << " " << clim.ppfd_max << " " << clim.ppfd << " " << clim.vpd << " " << clim.co2 << " " << clim.elv << " " << clim.swp << "\n"; 
}


void ClimateForcing::print_all(){
	for (int i = 0; i<t_met.size(); ++i){
		double t = t_met[i];
		int year = int(t);
		double month = (t-int(t))*12;
		std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
		std::cout << " | " << v_met[i].tc << " " << v_met[i].ppfd      << " " << v_met[i].vpd << "\n"; 
	}
}

double ClimateForcing::findMinDeltaT(std::vector <double> t_vector){
	int n = t_vector.size();
    std::vector <double> dts = diff(t_vector);
	auto delta = min_element(dts.begin(), dts.end()); 
    return *delta;
}


} // namespace env


