#include <iomanip>
#include <cmath>
#include <stdexcept>

#include "climate.h"
using namespace std;

namespace env{


int Climate::init(){
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
		// delta = t_met[n-1] - t_met[0] + (t_met[n-1]-t_met[0])/(n-1)*(1+1e-20); //+ 1e-20;
		delta = t_met[n-1] - t_met[0] + (t_met[1]-t_met[0]); //+ 1e-20;
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
			std::getline(fin_co2, line);
		
			std::stringstream lineStream(line);

			std::string cell;
			
			std::getline(lineStream, cell, ',');
			int year = as<int>(cell);

			std::getline(lineStream, cell, ',');
			double co2 = as<double>(cell);
			
			t_co2.push_back(year);
			v_co2.push_back(co2);
		}
	}
	
	return 0;

}


void Climate::set(double _tc, double _ppfd_max, double _ppfd, double _vpd, double _co2, double _elv, double _swp){
	clim.tc = _tc; 
	clim.ppfd_max = _ppfd_max; 
	clim.ppfd = _ppfd;
	clim.vpd = _vpd;
	clim.co2 = _co2;
	clim.elv = _elv;
	clim.swp = _swp;
};


Clim Climate::interp(Clim &clim_prev, Clim &clim_next){
	return clim_prev;
}

int Climate::id(double t){
	int id = (t - *t_met.begin())/delta*t_met.size();
	id = id % t_met.size();
	return id;
}


void Climate::updateClimate(double t){

	// double tadj = t;  // adjusted t to lie between limits of observed data
	// while(tadj < *t_met.begin()) tadj += delta;
		
	if (update_met){
		int n = t_met.size();
		if (n == 1){
			clim = v_met[0];
		}
		else{
			double tadj = t;  // adjusted t to lie between limits of observed data
			while(tadj < t_met[0])   tadj += delta;
			int idx_now = id(tadj);
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

int Climate::readNextLine_met(Clim &clim, double &t, std::ifstream& fin_met){

	std::string                line, cell;
	
	// READ line
	std::getline(fin_met, line);

	std::stringstream          lineStream(line);
	
	// PARSE line
	std::getline(lineStream, cell, ',');
	int year = as<int>(cell);
	
	std::getline(lineStream, cell, ',');
	int month = as<int>(cell);

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

	t = year + (month-1)/12.0;
		
	return 0;
}


void Climate::print(){
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

void Climate::print_line(double t){
	int year = int(t);
	double month = (t-int(t))*12;
	std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
	std::cout << " | T/Im/I/D/CO2/Z/Ps = " << clim.tc << " " << clim.ppfd_max << " " << clim.ppfd << " " << clim.vpd << " " << clim.co2 << " " << clim.elv << " " << clim.swp << "\n"; 
}


void Climate::print_all(){
	for (int i = 0; i<t_met.size(); ++i){
		double t = t_met[i];
		int year = int(t);
		double month = (t-int(t))*12;
		std::cout << "Climate at t = " << t << " (" << year << "/" << month << ")";
		std::cout << " | " << v_met[i].tc << " " << v_met[i].ppfd      << " " << v_met[i].vpd << "\n"; 
	}
}



} // namespace env


