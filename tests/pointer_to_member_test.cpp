// lexicographical_compare example
#include <iostream>     // std::cout, std::boolalpha
#include <algorithm>    // std::lexicographical_compare
#include <cctype>       // std::tolower
#include <vector>
#include <string>
#include <unordered_map>

class Traits{
	protected:
	static std::unordered_map<std::string, double Traits::*> members;
	
	public:
	double a;
	double b;
	
	Traits(double _a, double _b) : a(_a), b(_b){
	}
	
	double& operator[](const std::string& s){
		return this->*(members.find(s)->second);
	}
	
	// void set(const std::string& s, double value){
	//     this->*Traits::members.find(s)->second = value;
	// }
};

std::unordered_map<std::string, double Traits::*> Traits::members {
	{"a", &Traits::a},
	{"b", &Traits::b}
};        

int main () {

	Traits A(1,2);
	Traits B(3,4);
	
	std::cout << "A::a = " << A.a << " = " << A["a"] << '\n';
	std::cout << "A::b = " << A.b << " = " << A["b"] << '\n';
	std::cout << "B::a = " << B.a << " = " << B["a"] << '\n';
	std::cout << "B::b = " << B.b << " = " << B["b"] << '\n';

	if (A["a"] != 1) return 1;
	if (A["b"] != 2) return 1;
	if (B["a"] != 3) return 1;
	if (B["b"] != 4) return 1;

	A["a"] = 1.1;
	B["b"] = 4.4;

	std::cout << "A::a = " << A.a << " = " << A["a"] << '\n';
	std::cout << "A::b = " << A.b << " = " << A["b"] << '\n';
	std::cout << "B::a = " << B.a << " = " << B["a"] << '\n';
	std::cout << "B::b = " << B.b << " = " << B["b"] << '\n';

	if (A["a"] != 1.1) return 1;
	if (A["b"] != 2) return 1;
	if (B["a"] != 3) return 1;
	if (B["b"] != 4.4) return 1;


  return 0;
}
