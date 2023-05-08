// cppimport
#include <pybind11/pybind11.h>
#include "plantfate.h"

namespace py = pybind11;

PYBIND11_MODULE(simulator, m)
{
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<std::string>())
        .def("init", &Simulator::init)
        .def("simulate", &Simulator::simulate)
        .def("close", &Simulator::close);
}
/*

*/