// Example program
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

double positive_fmod(double x, double period){
    double res = fmod(x, period);
    if (res < 0) res += period;
    return res;
}


int main(){
{
  double t1 = 0, t2 = 1;
  std::cout << "fmod: " << t1 + fmod(0.2-t1, t2-t1) << "\n"; // expected 0.2
  std::cout << "fmod: " << t1 + fmod(0.8-t1, t2-t1) << "\n"; // expected 0.8
  std::cout << "fmod: " << t1 + fmod(-0.1 -t1, t2-t1) << "\n"; // want 0.9
  std::cout << "fmod: " << t1 + fmod(3.4 -t1, t2-t1) << "\n"; // want 0.4
  std::cout << "fmod: " << t1 + fmod(-3.2 -t1, t2-t1) << "\n"; // want 0.8

  std::cout << "remainder: " << t1 + remainder(0.2-t1, t2-t1) << "\n"; // expected 0.2
  std::cout << "remainder: " << t1 + remainder(0.8-t1, t2-t1) << "\n"; // expected 0.8
  std::cout << "remainder: " << t1 + remainder(-0.1 -t1, t2-t1) << "\n"; // want 0.9
  std::cout << "remainder: " << t1 + remainder(3.4 -t1, t2-t1) << "\n"; // want 0.4
  std::cout << "remainder: " << t1 + remainder(-3.2 -t1, t2-t1) << "\n"; // want 0.8
}
{
  double t1 = 0, t2 = 1;

  std::cout << "fmod+: " << t1 + positive_fmod(0.2-t1, t2-t1) << "\n"; // expected 0.2
  std::cout << "fmod+: " << t1 + positive_fmod(0.8-t1, t2-t1) << "\n"; // expected 0.8
  std::cout << "fmod+: " << t1 + positive_fmod(-0.1 -t1, t2-t1) << "\n"; // want 0.9
  std::cout << "fmod+: " << t1 + positive_fmod(3.4 -t1, t2-t1) << "\n"; // want 0.4
  std::cout << "fmod+: " << t1 + positive_fmod(-3.2 -t1, t2-t1) << "\n"; // want 0.8


  std::cout << "floor: " << t1 + 0.2-t1 - floor((0.2 -t1)/ (t2-t1)) << "\n"; // expected 0.2
  std::cout << "floor: " << t1 + 0.8-t1 - floor((0.8 -t1)/ (t2-t1)) << "\n"; // expected 0.8
  std::cout << "floor: " << t1 + -0.1-t1 - floor((-0.1 -t1)/ (t2-t1)) << "\n"; // want 0.9
  std::cout << "floor: " << t1 + 3.4-t1  - floor((3.4 -t1)/ (t2-t1)) << "\n"; // want 0.4
  std::cout << "floor: " << t1 + -3.2-t1 - floor((-3.2 -t1)/ (t2-t1)) << "\n"; // want 0.8
}
}
