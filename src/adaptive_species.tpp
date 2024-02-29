
template <class Model>
AdaptiveSpecies<Model>::AdaptiveSpecies(Model M, bool res) : Species<Model>(M) {
	int n = get_traits().size();
	fitness_gradient.resize(n, 0);
	trait_scalars.resize(n,1);
	trait_variance.resize(n, 0.01);
	for (int i=0; i<n; ++i) trait_names.push_back("T"+std::to_string(i));
	invasion_fitness = 0;
	r0 = 1;
	isResident = res;
	r0_hist.set_interval(200);
}


template <class Model>
void AdaptiveSpecies<Model>::set_traits(std::vector<double> tvec){
	this->boundaryCohort.set_evolvableTraits(tvec);
	for (auto& c : this->cohorts) c.set_evolvableTraits(tvec);
}


template <class Model>
std::vector<double> AdaptiveSpecies<Model>::get_traits(){
	return this->boundaryCohort.get_evolvableTraits();
}


template <class Model>
void AdaptiveSpecies<Model>::createVariants(Model M){
	std::vector<double> traits = get_traits();
	for (int i=0; i<traits.size(); ++i){
		std::vector<double> traits_mutant = traits;
		traits_mutant[i] += fg_dx * trait_scalars[i];
		AdaptiveSpecies<Model> * m = new AdaptiveSpecies<Model>(*this); // copy-construct the variant, then edit traits
		m->set_traits(traits_mutant);
		m->isResident = false;
		m->probes.clear(); // mutants dont have probes
		m->species_name += "_probe" + std::to_string(i);
		probes.push_back(m);
	}
}


template <class Model>
void AdaptiveSpecies<Model>::evolveTraits(double dt){
	if (!isResident) return;

	std::vector<double> traits_res = get_traits();
	std::vector<double> dx(traits_res.size());

	for (int i=0; i<dx.size(); ++i){
		double dx_norm = fitness_gradient[i]*trait_variance[i]*dt; // Landes equation for normalized trait
		dx[i] = trait_scalars[i] * dx_norm;
	}

	for (int i=0; i<dx.size(); ++i) traits_res[i] += dx[i];
	set_traits(traits_res);

	for (auto m : probes){
		std::vector<double> traits_mut = m->get_traits();
		for (int j=0; j<traits_mut.size(); ++j) traits_mut[j] += dx[j];
		m->set_traits(traits_mut);
	}
}


template <class Model>
void AdaptiveSpecies<Model>::calcFitnessGradient(){
	if (!isResident) return;

	fitness_gradient.clear();
	for (auto m : probes){
		double grad = (m->r0_hist.get()-r0_hist.get()) / fg_dx;
		fitness_gradient.push_back(grad);
		//cout << "   " << m->invasion_fitness << " " << spp->invasion_fitness << " " << *spp->fitness_gradient.rbegin() << "\n";
	}
}


template <class Model>
void AdaptiveSpecies<Model>::print_extra(){
	std::cout << "Name: " << species_name << "\n";
	std::cout << "Resident: " << ((isResident)? "Yes" : "No") << "\n";
	std::cout << "Probes: ";
	for (auto x : probes) std::cout << x->species_name << " ";
	std::cout << "\n";
	std::cout << "Fitness gradient: ";
	for (auto x : fitness_gradient) std::cout << x << " ";
	std::cout << "\n";
	std::cout << "Trait variance: ";
	for (auto x : trait_variance) std::cout << x << " ";
	std::cout << "\n";
	std::cout << "Trait scalars: ";
	for (auto x : trait_scalars) std::cout << x << " ";
	std::cout << "\n";     // these scalars will be applied to fg_dx
}

// Changelog:
// v2: Save plant parameters also from boundary cohort, so that ini file is not needed during restore (for plant parameters)
template <class Model>
void AdaptiveSpecies<Model>::save(std::ostream &fout){
	fout << "AdaptiveSpecies<T>::v2\n";
	
	// save species-level data
	fout << std::make_tuple(
		      fg_dx
			, std::quoted(species_name)
			, isResident
			, t_introduction
			, invasion_fitness
			, r0);
	fout << '\n';

	fout << fitness_gradient
	     << trait_variance
		 << trait_scalars
		 << trait_names;  // trait names dont contain whitespaces, so safe to write this way
	
	// MovingAveragers
	seeds_hist.save(fout);
	r0_hist.save(fout);

	// save both par and traits from boundary cohort
	auto& C = this->getCohort(-1);
	C.par.save(fout); 
	C.traits.save(fout); 

	Species<Model>::save(fout);
}


template <class Model>
void AdaptiveSpecies<Model>::restore(std::istream &fin){
	std::cout << "Restoring AdaptiveSpecies<Model>...\n";
	std::string s; fin >> s; // discard version number
	assert(s == "AdaptiveSpecies<T>::v2");

	if (configfile_for_restore == "") throw std::runtime_error("Config file has not been set");

	// restore species-level data
	fin >> fg_dx
		>> std::quoted(species_name)
		>> isResident
		>> t_introduction
		>> invasion_fitness
		>> r0;
	
	fin >> fitness_gradient
	    >> trait_variance
		>> trait_scalars
		>> trait_names;

	// moving averagers go here
	seeds_hist.restore(fin);
	r0_hist.restore(fin);

	// Create a Model object and restore all individual properties to this object
	// This will be used to copy-construct the species
	auto& C = this->getCohort(-1);
	//C.initFromFile(configfile_for_restore);  
	C.par.restore(fin);
	C.traits.restore(fin);
	C.init(C.par, C.traits);

	C.traits.save(std::cout); std::cout.flush(); // For debug only

	Species<Model>::restore(fin);
}

