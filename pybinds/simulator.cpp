// cppimport
#include <pybind11/pybind11.h>
#include "plantfate.h"
#include "climate.h"
#include "pspm_interface.h"

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}


PYBIND11_MODULE(simulator, m)
{
    m.def("add", &add);
    py::class_<env::Clim>(m, "Clim", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("co2", &env::Clim::co2);

    py::class_<Simulator>(m, "Simulator", py::dynamic_attr())
        .def(py::init<std::string>())
        .def_readwrite("paramsFile", &Simulator::paramsFile)
		.def_readwrite("parent_dir", &Simulator::parent_dir)
		.def_readwrite("expt_dir", &Simulator::expt_dir)
		.def_readwrite("met_file", &Simulator::met_file)
		.def_readwrite("co2_file", &Simulator::co2_file)
		.def_readwrite("save_state", &Simulator::save_state)
		.def_readwrite("state_outfile", &Simulator::state_outfile)
		.def_readwrite("config_outfile", &Simulator::config_outfile)
		.def_readwrite("continueFrom_stateFile", &Simulator::continueFrom_stateFile)
		.def_readwrite("continueFrom_configFile", &Simulator::continueFrom_configFile)
		.def_readwrite("continuePrevious", &Simulator::continuePrevious)
		.def_readwrite("evolve_traits", &Simulator::evolve_traits)
		.def_readwrite("y0", &Simulator::y0)
		.def_readwrite("yf", &Simulator::yf)
		.def_readwrite("ye", &Simulator::ye);

//     py::class_<PSPM_Plant, plant::Plant>(m, "PSPM_Plant")
//         .def(py::init<>())
//         .def_readwrite("t_birth", &PSPM_Plant::t_birth)
//         .def_readwrite("varnames", &PSPM_Plant::varnames)
//         .def_readwrite("statevarnames", &PSPM_Plant::statevarnames)
//         .def_readwrite("nrc", &PSPM_Plant::nrc)
//         .def_readwrite("ndc", &PSPM_Plant::ndc)
//         .def_readwrite("nbc", &PSPM_Plant::nbc)
//         .def("set_size", &PSPM_Plant::set_size)
//         .def("init_density", &PSPM_Plant::init_density)
//         .def("preCompute", &PSPM_Plant::preCompute)
//         .def("afterStep", &PSPM_Plant::afterStep)
//         .def("establishmentProbability", &PSPM_Plant::establishmentProbability)
//         .def("growthRate", &PSPM_Plant::growthRate)
//         .def("mortalityRate", &PSPM_Plant::mortalityRate)
//         .def("birthRate", &PSPM_Plant::birthRate)
//         .def("init_state", &PSPM_Plant::init_state)
//         .def("set_state", &PSPM_Plant::set_state)
//         .def("get_state", &PSPM_Plant::get_state)
//         .def("get_rates", &PSPM_Plant::get_rates)
//         .def("print", &PSPM_Plant::print)
//         .def("save", &PSPM_Plant::save)
//         .def("restore", &PSPM_Plant::restore);
}



// PYBIND11_MODULE(simulator, m)
// {
    
//         // .def("simulate", &Simulator::simulate)
//         // .def("close", &Simulator::close);
// // }
// 
/*
<%
setup_pybind11(cfg)
cfg['include_dirs'] = ['../src', '../inst/include', '/opt/homebrew/Cellar/pybind11/2.10.3/include', '../../phydro/inst/include', '../../libpspm/include']
cfg['libraries'] = ['../../libpspm/lib']
cfg['extra_compile_args'] = ['-O3', '-g', '-fPIC', '-pg', '-std=c++17', '-Wall', '-Wextra', '-DPHYDRO_ANALYTICAL_ONLY']
%>
*/