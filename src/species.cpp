#include <species.h>


Species_Base::~Species_Base(){
}

void Species_Base::print_extra(){
}
	
// FIXME JJ: I am removing this because `cohort` is not the correct generic name. 
//           We use it internally to define eithe cohort/gridcell/etc, 
//           but the user need not be concerned with this 
// int Species_Base::cohortsize(){
// 	return J;
// }

//template<class Model>
//int Species<Model>::addVar(std::string name, int stride, int offset){
//    varnames.push_back(name);
//     strides.push_back(stride);
//    offsets.push_back(offset);
//}


//template<class Model>
//void Species<Model>::clearVars(){
//    varnames.clear();
//    strides.clear();
//    offsets.clear();
//}


//template<class Model>
//void Species<Model>::set_model(Model *M){
//    mod = M;
//}


void Species_Base::set_inputBirthFlux(double b){
	birth_flux_in = b;
}


//template<class Model>
//double Species<Model>::set_iStateVariables(std::vector<std::string> names){
//    varnames_extra = names;
//}

void Species_Base::set_bfin_is_u0in(bool flag){
	bfin_is_u0in = flag;
}


int Species_Base::xsize(){
	return J;
}


int Species_Base::stateSizeTotal(){
	return (istate_size + 1 + n_accumulators)*J;
	//      ^ x           ^ u      ^ accumulators
}


// void Species_Base::save(std::ofstream &fout){
// 	// fout << "Species_Base::v1\n";

// 	// fout << std::make_tuple(
// 	// 	J
// 	//   , n_extra_statevars
// 	//   , noff_abm
// 	//   , birth_flux_in
// 	//   , bfin_is_u0in
// 	//   , xb);
// 	// fout << '\n';
// 	// fout << X << x << h;
// }

// void Species_Base::restore(std::ifstream &fin){
// 	// std::cout << "Restoring Species_Base" << "\n";
// 	// std::string s; fin >> s; // version number (discard)
// 	// assert(s == "Species_Base::v1");
// 	// fin >> J	
// 	//     >> n_extra_statevars
// 	//     >> noff_abm
// 	//     >> birth_flux_in
// 	//     >> bfin_is_u0in
// 	//     >> xb;
	
// 	// fin >> X >> x >> h;
// }



