#include <pybind11/pybind11.h>
#include <environment_base.h>
#include <solver.h>

class PyEnvironmentBase : public EnvironmentBase {
public:
    /* Inherit the constructors */
    using EnvironmentBase::EnvironmentBase;

	~PyEnvironmentBase() override = default;

    /* Trampoline (need one for each virtual function) */
    void computeEnv(double t, Solver * sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt) override {
        PYBIND11_OVERRIDE_PURE(
            void, /* Return type */
            EnvironmentBase,      /* Parent class */
            computeEnv,          /* Name of function in C++ (must match Python name) */
            t, sol, S, dSdt      /* Argument(s) */
        );
    }
};
