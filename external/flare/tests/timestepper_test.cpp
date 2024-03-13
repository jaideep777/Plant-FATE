#include "time_stepper.h"
using namespace std;

int main(){

	flare::TimeStepper stepper;

	stepper.set_units("months since 2000-01-15 0:0:0");

	double j1 = stepper.to_julian(1);
	cout << "1 month since 2000-01-15 0:0:0 is: " << flare::julian_to_datestring(j1) << '\n';

	stepper.set_units("days since 2000-01-15 0:0:0");

	double j2 = stepper.to_julian(31);
	cout << "31 days since 2000-01-15 0:0:0 is: " << flare::julian_to_datestring(j2) << '\n';

	double t;
	stepper.j_current = flare::datestring_to_julian("2021-03-15");
	t = stepper.step(31);
	cout << "2021-03-15 0:0:0 + 31 days is: " << flare::julian_to_datestring(t) << '\n';
	cout << "stepper time in my units: " << stepper.to_stepperUnits(stepper.j_current) << '\n';
	if (fabs(stepper.to_stepperUnits(stepper.j_current) - 7761) > 1e-6) return 1;

	t = stepper.step(30);
	cout << "2021-03-15 0:0:0 + 31 + 30 days is: " << flare::julian_to_datestring(t) << '\n';
	cout << "stepper time in my units: " << stepper.to_stepperUnits(stepper.j_current) << '\n';
	if (fabs(stepper.to_stepperUnits(stepper.j_current) - 7791) > 1e-6) return 1;

	stepper.j_current = stepper.to_julian(1); // 1 day since 2000-01-15 = 2000-01-16
	t = stepper.step(5); // 5 more days since then = 2000-01-21
	cout << "new time (expected 2000-01-21): " << flare::julian_to_datestring(t) << '\n';
	cout << "stepper time in my units: " << stepper.to_stepperUnits(stepper.j_current) << '\n';
	if (fabs(stepper.to_stepperUnits(stepper.j_current) - 6) > 1e-6) return 1;

	return 0;
}
