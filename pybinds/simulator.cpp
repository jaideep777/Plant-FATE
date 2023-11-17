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
#include "climate_forcing.cpp"
#include "climate_input.cpp"

namespace py = pybind11;

PYBIND11_MODULE(plantFATE, m)
{
    py::class_<env::Clim>(m, "Clim", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("tc", &env::Clim::tc)
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

	// py::class_<env::ClimateForcing>(m, "ClimateForcing")
	// 	.def(py::init<>())
	// 	.def("init", &env::ClimateForcing::init)
	// 	.def("set", &env::ClimateForcing::set)
	// 	.def("update_tc", &env::ClimateForcing::update_tc)
	// 	.def("update_ppfd_max", &env::ClimateForcing::update_ppfd_max)
	// 	.def("update_ppfd", &env::ClimateForcing::update_ppfd)
	// 	.def("update_vpd", &env::ClimateForcing::update_vpd)
	// 	.def("update_elv", &env::ClimateForcing::update_elv)
	// 	.def("update_swp", &env::ClimateForcing::update_swp)
	// 	.def("id", &env::ClimateForcing::id)
	// 	.def("idx_of", &env::ClimateForcing::idx_of)
	// 	.def("updateClimate", &env::ClimateForcing::updateClimate)
	// 	.def("print", &env::ClimateForcing::print)
	// 	.def_readwrite("clim", &env::ClimateForcing::clim)
	// 	.def_readwrite("t_now", &env::ClimateForcing::t_now)
	// 	.def_readwrite("t_met", &env::ClimateForcing::t_met)
	// 	.def_readwrite("v_met", &env::ClimateForcing::v_met)
	// 	.def_readwrite("interpolate", &env::ClimateForcing::interpolate)
	// 	.def_readwrite("update_met", &env::ClimateForcing::update_met)
	// 	.def_readwrite("update_co2", &env::ClimateForcing::update_co2);
	

	py::class_<env::ClimateInput>(m, "ClimateInput")
		.def(py::init<>())
		.def(py::init<env::Clim&, double, double>())
		.def("updateEnvironment", &env::ClimateInput::updateEnvironment)
		.def("updateClim", &env::ClimateInput::updateClim)
		.def("print_line", &env::ClimateInput::print_line)
		.def_readwrite("currentClim", &env::ClimateInput::currentClim)
		.def_readwrite("weightedAveClim", &env::ClimateInput::weightedAveClim)
		.def_readwrite("tcurrent", &env::ClimateInput::tcurrent)
		.def_readwrite("ave_window", &env::ClimateInput::ave_window);

	py::class_<PSPM_Dynamic_Environment, env::LightEnvironment, env::ClimateInput>(m, "PSPM_Dynamic_Environment")
		.def(py::init<>());

	py::class_<ErgodicEnvironment, env::LightEnvironment, env::ClimateInput>(m, "ErgodicEnvironment")
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

	py::class_<Simulator>(m, "Simulator", py::dynamic_attr())
		.def(py::init<std::string>())
		.def("init", py::overload_cast<double, double>(&Simulator::init), py::arg("tstart"), py::arg("tend")) 
		.def("init", py::overload_cast<double>(&Simulator::init), py::arg("tstart")) 
		.def("init", py::overload_cast<double, env::Clim&>(&Simulator::init), py::arg("tstart"), py::arg("initClimObj")) 
		// .def("init", &Simulator::init, (void(Simulator::*)(double, double))&Simulator::init, py::arg("tstart", "tend"))
		// .def("init", &Simulator::init, (void(Simulator::*)(double))&Simulator::init, py::arg("tstart"))
		// .def("init", &Simulator::init, (void(Simulator::*)(double, &env::Climate))&Simulator::init, py::arg("tstart", "climobj"))
		.def("simulate", &Simulator::simulate)
		.def("simulate_to", &Simulator::simulate_to)
		.def("simulate_step", &Simulator::simulate_step)
		.def("close", &Simulator::close)
		// .def("set_metFile", &Simulator::set_metFile)
		// .def("set_co2File", &Simulator::set_co2File)
		.def("update_environment", &Simulator::update_environment)
		.def_readwrite("E", &Simulator::E)
        .def_readwrite("paramsFile", &Simulator::paramsFile)
		.def_readwrite("parent_dir", &Simulator::parent_dir)
		.def_readwrite("expt_dir", &Simulator::expt_dir)
		.def_readwrite("save_state", &Simulator::save_state)
		.def_readwrite("state_outfile", &Simulator::state_outfile)
		.def_readwrite("config_outfile", &Simulator::config_outfile)
		.def_readwrite("continueFrom_stateFile", &Simulator::continueFrom_stateFile)
		.def_readwrite("continueFrom_configFile", &Simulator::continueFrom_configFile)
		.def_readwrite("continuePrevious", &Simulator::continuePrevious)
		.def_readwrite("evolve_traits", &Simulator::evolve_traits)
		.def_readwrite("y0", &Simulator::y0)
		.def_readwrite("yf", &Simulator::yf)
		.def_readwrite("ye", &Simulator::ye)
		.def_readwrite("props", &Simulator::props)
		.def_readwrite("cwm", &Simulator::cwm)
		.def_readwrite("tcurrent", &Simulator::tcurrent)
		.def_readwrite("delta_T", &Simulator::delta_T)
		.def_readwrite("T_invasion", &Simulator::T_invasion);
		// .def_readwrite("traits0", &Simulator::traits0)
		// .def_readwrite("par0", &Simulator::par0);

	py::class_<plant::PlantTraits>(m, "PlantTraits")
		.def(py::init<>())
		.def("print", &plant::PlantTraits::print)
		// .def_readwrite("lma", &plant::PlantTraits::lma)
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

py::class_<EmergentProps>(m, "EmergentProps")
		.def(py::init<>())
		.def_readwrite("gpp", &EmergentProps::gpp)
		.def_readwrite("npp", &EmergentProps::npp)
		.def_readwrite("resp_auto", &EmergentProps::resp_auto)
		.def_readwrite("trans", &EmergentProps::trans)
		.def_readwrite("gs", &EmergentProps::gs)
		.def_readwrite("lai", &EmergentProps::lai)
		.def_readwrite("leaf_mass", &EmergentProps::leaf_mass)
		.def_readwrite("stem_mass", &EmergentProps::stem_mass)
		.def_readwrite("croot_mass", &EmergentProps::croot_mass)
		.def_readwrite("froot_mass", &EmergentProps::froot_mass)
		.def_readwrite("cc_est", &EmergentProps::cc_est)
		.def_readwrite("lai_vert", &EmergentProps::lai_vert);

py::class_<SpeciesProps>(m, "SpeciesProps")
		.def(py::init<>())
		.def_readwrite("n_ind", &SpeciesProps::n_ind)
		.def_readwrite("biomass", &SpeciesProps::biomass)
		.def_readwrite("ba", &SpeciesProps::ba)
		.def_readwrite("canopy_area", &SpeciesProps::canopy_area)
		.def_readwrite("height", &SpeciesProps::height)
		.def_readwrite("lma", &SpeciesProps::lma)
		.def_readwrite("p50", &SpeciesProps::p50)
		.def_readwrite("hmat", &SpeciesProps::hmat)
		.def_readwrite("wd", &SpeciesProps::wd)
		.def_readwrite("gs", &SpeciesProps::gs)
		.def_readwrite("vcmax", &SpeciesProps::vcmax)
		.def_readwrite("n_ind_vec", &SpeciesProps::n_ind_vec)
		.def_readwrite("biomass_vec", &SpeciesProps::biomass_vec)
		.def_readwrite("ba_vec", &SpeciesProps::ba_vec)
		.def_readwrite("canopy_area_vec", &SpeciesProps::canopy_area_vec)
		.def_readwrite("height_vec", &SpeciesProps::height_vec)
		.def_readwrite("vcmax_vec", &SpeciesProps::vcmax_vec)
		.def_readwrite("lma_vec", &SpeciesProps::lma_vec)
		.def_readwrite("p50_vec", &SpeciesProps::p50_vec)
		.def_readwrite("hmat_vec", &SpeciesProps::hmat_vec)
		.def_readwrite("wf_vec", &SpeciesProps::wd_vec);



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