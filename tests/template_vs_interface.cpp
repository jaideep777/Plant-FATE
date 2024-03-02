#include <vector>
#include <iostream>
using namespace std;


class Species{
	public:
	virtual double g(){
		return 0;
	}
};

class IndianSpecies : public Species{
	virtual double g(){
		return 1;
	}
};


template<class T>
class TemplateSolver{
	public:
	vector<T*> vec;

	void print(){
		for (int i=0; i<vec.size(); ++i){
			cout << vec[i]->g() << " "; 
		}
		cout << endl;
	}

};

class InterfaceSolver{
	public:
	vector<Species*> vec;
	
	void print(){
		for (int i=0; i<vec.size(); ++i){
			cout << vec[i]->g() << " "; 
		}
		cout << endl;
	}

};

int main(){
	Species s;
	IndianSpecies i;

	InterfaceSolver Si;
	Si.vec.push_back(&s);
	Si.vec.push_back(&i);
	Si.print();

	TemplateSolver<Species> St;
	St.vec.push_back(&s);
	St.vec.push_back(&i);
	St.print();



}
