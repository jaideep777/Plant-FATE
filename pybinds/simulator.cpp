// cppimport
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "pybind_override.cpp"
#include "environment_base.h"
#include "plantfate.cpp"
#include "light_environment.cpp"
#include "climate.cpp"
#include "pspm_interface.cpp"
#include "pspm_dynamic_environment.cpp"
#include "plant_geometry.cpp"
#include "assimilation.cpp"
#include "plant.cpp"
#include "state_restore.cpp"
#include "trait_reader.h"
#include "treelife.cpp" 
#include "trait_evolution.h"
#include "community_properties.cpp"
#include "plant_params.h"

namespace py = pybind11;

PYBIND11_MODULE(plantFATE, m)
{
    py::class_<env::Clim>(m, "Clim", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("co2", &env::Clim::tc)
		.def_readwrite("ppfd_max", &env::Clim::ppfd_max)
		.def_readwrite("ppfd", &env::Clim::ppfd)
		.def_readwrite("vpd", &env::Clim::vpd)
		.def_readwrite("co2", &env::Clim::co2)
		.def_readwrite("elv", &env::Clim::elv)
		.def_readwrite("swp", &env::Clim::swp);

	py::class_<EnvironmentBase, PyEnvironmentBase>(m, "EnvironmentBase")
	    .def(py::init<>())
		.def("computeEnv", &EnvironmentBase::computeEnv);

	py::class_<env::LightEnvironment, EnvironmentBase>(m, "LightEnvironment")
	    .def(py::init<>())
		.def("computeEnv", &env::LightEnvironment::computeEnv)
		.def("print", &env::LightEnvironment::print)
		.def_readwrite("use_ppa", &env::LightEnvironment::use_ppa)
		.def_readwrite("n_layers", &env::LightEnvironment::n_layers)
		.def_readwrite("z_star", &env::LightEnvironment::z_star)
		.def_readwrite("fapar_tot", &env::LightEnvironment::fapar_tot)
		.def_readwrite("canopy_openness", &env::LightEnvironment::canopy_openness);

	py::class_<env::Climate>(m, "Climate")
		.def(py::init<>())
		.def("print", &env::Climate::print)
		.def_readwrite("clim", &env::Climate::clim);
	
	py::class_<PSPM_Dynamic_Environment, env::LightEnvironment, env::Climate>(m, "PSPM_Dynamic_Environment")
		.def(py::init<>());

	py::class_<ErgodicEnvironment, env::LightEnvironment, env::Climate>(m, "ErgodicEnvironment")
		.def(py::init<>())
		.def("print", &ErgodicEnvironment::print);

	py::class_<LifeHistoryOptimizer>(m,"LifeHistoryOptimizer")
		.def(py::init<std::string>())
		// .def("set_metFile", &LifeHistoryOptimizer::set_metFile)
		// .def("set_co2File", &LifeHistoryOptimizer::set_co2File)
		// .def("get_header", &LifeHistoryOptimizer::get_header)
		// .def("get_state", &LifeHistoryOptimizer::get_state)
		.def("set_traits", &LifeHistoryOptimizer::set_traits)
		.def("get_traits", &LifeHistoryOptimizer::get_traits)
		.def("init", &LifeHistoryOptimizer::init)
		// .def("printMeta", &LifeHistoryOptimizer::printMeta)
		.def("calcFitness", &LifeHistoryOptimizer::calcFitness)
		.def("grow_for_dt", &LifeHistoryOptimizer::grow_for_dt);
		// .def_readwrite("params_file", &LifeHistoryOptimizer::params_file);
		// .def_readwrite("env", &LifeHistoryOptimizer::C);
		// .def_readwrite("traits0", &LifeHistoryOptimizer::traits0)
		// .def_readwrite("par0", &LifeHistoryOptimizer::par0);

	py::class_<Patch>(m, "Patch")
		.def(py::init<std::string>())
		.def("init", &Patch::init)
		.def("simulate", &Patch::simulate)
		.def("close", &Patch::close)
		// .def("set_metFile", &Patch::set_metFile)
		// .def("set_co2File", &Patch::set_co2File)
		// .def_readwrite("E", &Patch::E)
        .def_readwrite("paramsFile", &Patch::paramsFile)
		.def_readwrite("parent_dir", &Patch::parent_dir)
		.def_readwrite("expt_dir", &Patch::expt_dir)
		// .def_readwrite("met_file", &Patch::met_file)
		// .def_readwrite("co2_file", &Patch::co2_file)
		.def_readwrite("save_state", &Patch::save_state)
		.def_readwrite("state_outfile", &Patch::state_outfile)
		.def_readwrite("config_outfile", &Patch::config_outfile)
		.def_readwrite("continueFrom_stateFile", &Patch::continueFrom_stateFile)
		.def_readwrite("continueFrom_configFile", &Patch::continueFrom_configFile)
		.def_readwrite("continuePrevious", &Patch::continuePrevious)
		.def_readwrite("evolve_traits", &Patch::evolve_traits)
		.def_readwrite("y0", &Patch::y0)
		.def_readwrite("yf", &Patch::yf)
		.def_readwrite("ye", &Patch::ye);
		// .def_readwrite("traits0", &Patch::traits0)
		// .def_readwrite("par0", &Patch::par0);

	py::class_<plant::PlantTraits>(m, "PlantTraits")
		.def(py::init<>())
		.def("print", &plant::PlantTraits::print)
		.def_readwrite("lma", &plant::PlantTraits::lma)
		.def_readwrite("zeta", &plant::PlantTraits::zeta)
		.def_readwrite("fcr", &plant::PlantTraits::fcr)
		.def_readwrite("hmat", &plant::PlantTraits::hmat)
		.def_readwrite("fhmat", &plant::PlantTraits::fhmat)
		.def_readwrite("seed_mass", &plant::PlantTraits::seed_mass)
		.def_readwrite("wood_density", &plant::PlantTraits::wood_density)
		.def_readwrite("p50_xylem", &plant::PlantTraits::p50_xylem)
		.def_readwrite("K_leaf", &plant::PlantTraits::K_leaf)
		.def_readwrite("K_xylem", &plant::PlantTraits::K_xylem)
		.def_readwrite("b_leaf", &plant::PlantTraits::b_leaf)
		.def_readwrite("b_xylem", &plant::PlantTraits::b_xylem)
		// .def_readwrite("ll", &plant::PlantTraits::ll)
		.def_readwrite("p50_leaf", &plant::PlantTraits::p50_leaf);

	py::class_<plant::PlantParameters>(m, "PlantParameters")
		.def(py::init<>())
		.def("print", &plant::PlantParameters::print)
		.def_readwrite("cD0", &plant::PlantParameters::cD0)
		// .def_readwrite("eD0", &plant::PlantParameters::eD0)
		.def_readwrite("cD1", &plant::PlantParameters::cD1)
		.def_readwrite("m_alpha", &plant::PlantParameters::m_alpha)
		.def_readwrite("m_beta", &plant::PlantParameters::m_beta)
		.def_readwrite("m_gamma", &plant::PlantParameters::m_gamma);
		// .def_readwrite("eWD_alpha", &plant::PlantParameters::eWD_alpha)
		// .def_readwrite("eWD_gamma", &plant::PlantParameters::eWD_gamma)
		// .def_readwrite("cDW0", &plant::PlantParameters::cDW0)
		// .def_readwrite("eWD", &plant::PlantParameters::eWD);

	

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
    
//         // .def("simulate", &Patch::simulate)
//         // .def("close", &Patch::close);
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