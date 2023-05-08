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
<%
setup_pybind11(cfg)
cfg['include_dirs'] = ['../src', '../inst/include', '/opt/homebrew/Cellar/pybind11/2.10.3/include', '../../phydro/inst/include', '../../libpspm/include']
cfg['libraries'] = ['../../libpspm/lib']
cfg['extra_compile_args'] = ['-O3', '-g', '-fPIC', '-pg', '-std=c++17', '-Wall', '-Wextra', '-DPHYDRO_ANALYTICAL_ONLY']
%>
*/