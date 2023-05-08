#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pspm_interface.h"

namespace py = pybind11;

PYBIND11_MODULE(pspm_interface, m)
{
    py::class_<PSPM_Plant, plant::Plant>(m, "PSPM_Plant")
        .def(py::init<>())
        .def_readwrite("t_birth", &PSPM_Plant::t_birth)
        .def_readwrite("varnames", &PSPM_Plant::varnames)
        .def_readwrite("statevarnames", &PSPM_Plant::statevarnames)
        .def_readwrite("nrc", &PSPM_Plant::nrc)
        .def_readwrite("ndc", &PSPM_Plant::ndc)
        .def_readwrite("nbc", &PSPM_Plant::nbc)
        .def("set_size", &PSPM_Plant::set_size)
        .def("init_density", &PSPM_Plant::init_density)
        .def("preCompute", &PSPM_Plant::preCompute)
        .def("afterStep", &PSPM_Plant::afterStep)
        .def("establishmentProbability", &PSPM_Plant::establishmentProbability)
        .def("growthRate", &PSPM_Plant::growthRate)
        .def("mortalityRate", &PSPM_Plant::mortalityRate)
        .def("birthRate", &PSPM_Plant::birthRate)
        .def("init_state", &PSPM_Plant::init_state)
        .def("set_state", &PSPM_Plant::set_state)
        .def("get_state", &PSPM_Plant::get_state)
        .def("get_rates", &PSPM_Plant::get_rates)
        .def("print", &PSPM_Plant::print)
        .def("save", &PSPM_Plant::save)
        .def("restore", &PSPM_Plant::restore);
}
