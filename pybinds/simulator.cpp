// cppimport
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "pybind_override.cpp"
#include "environment_base.h"
#include "plantfate.h"
#include "light_environment.cpp"
#include "climate.h"
#include "pspm_interface.h"
#include "test_environment.cpp"

namespace py = pybind11;


// class TestEnvironment: public EnvironmentBase {
// public:

//     bool use_ppa = false;

//     //PPA Environment:
// 	int n_layers;
// 	double total_crown_area;
// 	std::vector<double> z_star;
// 	std::vector<double> fapar_tot;
// 	std::vector<double> canopy_openness;


//     TestEnvironment(){
// 		n_layers = 0;
// 		z_star = {0};
// 		fapar_tot = {0};
// 		canopy_openness = {1};
// 		z_star.reserve(20);
// 	    canopy_openness.reserve(20);
// 	}

//     void computeEnv(double t, Solver * sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt) {
// 		std::cout<< "working" << std::endl;
// 	}
	
// 	void print(){
// 		std::cout << "PPA:" << "\n";
// 		std::cout << "z* (" << n_layers << ") = "; for (auto z:z_star) std::cout << z << " "; 
// 		std::cout << "\nfapar layer* = "; for (auto z:fapar_tot) std::cout << z << " "; 
// 		std::cout << "\ncanopy openness* = "; for (auto z:canopy_openness) std::cout << z << " "; 
// 		std::cout << "\n";
// 	}

// };


PYBIND11_MODULE(simulator, m)
{
    py::class_<env::Clim>(m, "Clim", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("co2", &env::Clim::co2);

	py::class_<EnvironmentBase, PyEnvironmentBase>(m, "EnvironmentBase")
	    .def(py::init<>())
		.def("computeEnv", &EnvironmentBase::computeEnv);

    py::class_<TestEnvironment, EnvironmentBase>(m, "TestEnvironment")
	    .def(py::init<>())
		.def("computeEnv", &TestEnvironment::computeEnv)
		.def("print", &TestEnvironment::print)
		.def_readwrite("use_ppa", &TestEnvironment::use_ppa)
		.def_readwrite("n_layers", &TestEnvironment::n_layers)
		.def_readwrite("z_star", &TestEnvironment::z_star)
		.def_readwrite("fapar_tot", &TestEnvironment::fapar_tot)
		.def_readwrite("canopy_openness", &TestEnvironment::canopy_openness);


	// py::class_<LightEnvironment, EnvironmentBase>(m, "LightEnvironment")
	//     .def(py::init<>())
	// 	.def("computeEnv", &LightEnvironment::computeEnv)
	// 	.def("print", &LightEnvironment::print)
	// 	.def_readwrite("use_ppa", &LightEnvironment::use_ppa)
	// 	.def_readwrite("n_layers", &LightEnvironment::n_layers)
	// 	.def_readwrite("z_star", &LightEnvironment::z_star)
	// 	.def_readwrite("fapar_tot", &LightEnvironment::fapar_tot)
	// 	.def_readwrite("canopy_openness", &LightEnvironment::canopy_openness);

		
	py::class_<env::Climate> Climate(m, "Climate");
	// py::class_<PSPM_Dynamic_Environment, EnvironmentBase, env::LightEnvironment, env::Climate>(m, "PSPM_Dynamic_Environment");
	// py::class_<PSPM_Dynamic_Environment>(m, "PSPM_Dynamic_Environment");
    // py::class_<Simulator> Simulator(m, "Simulator");
        // .def(py::init<std::string>())
        // .def_readwrite("paramsFile", &Simulator::paramsFile)
		// .def_readwrite("parent_dir", &Simulator::parent_dir)
		// .def_readwrite("expt_dir", &Simulator::expt_dir)
		// .def_readwrite("met_file", &Simulator::met_file)
		// .def_readwrite("co2_file", &Simulator::co2_file)
		// .def_readwrite("save_state", &Simulator::save_state)
		// .def_readwrite("state_outfile", &Simulator::state_outfile)
		// .def_readwrite("config_outfile", &Simulator::config_outfile)
		// .def_readwrite("continueFrom_stateFile", &Simulator::continueFrom_stateFile)
		// .def_readwrite("continueFrom_configFile", &Simulator::continueFrom_configFile)
		// .def_readwrite("continuePrevious", &Simulator::continuePrevious)
		// .def_readwrite("evolve_traits", &Simulator::evolve_traits)
		// .def_readwrite("y0", &Simulator::y0)
		// .def_readwrite("yf", &Simulator::yf)
		// .def_readwrite("ye", &Simulator::ye);

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
};



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