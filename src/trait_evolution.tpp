
template <class Model>
MySpecies<Model>::MySpecies(Model M, bool res) : Species<Model>(M) {
	int n = get_traits().size();
	fitness_gradient.resize(n, 0);
	trait_variance.resize(n, 0.01);
	invasion_fitness = 0;
	r0 = 1;
	isResident = res;
}

template <class Model>
void MySpecies<Model>::set_traits(std::vector<double> tvec){
	this->boundaryCohort.set_traits(tvec);
	for (auto& c : this->cohorts) c.set_traits(tvec);
}

template <class Model>
std::vector<double> MySpecies<Model>::get_traits(){
	return this->boundaryCohort.get_traits();
}


template <class Model>
void MySpecies<Model>::createVariants(Model M){
	std::vector<double> traits = get_traits();
	for (int i=0; i<traits.size(); ++i){
		std::vector<double> traits_mutant = traits;
		traits_mutant[i] += fg_dx * trait_scalars[i];
		MySpecies<Model> * m = new MySpecies<Model>(*this); // copy-construct the variant, then edit traits
		m->set_traits(traits_mutant);
		m->isResident = false;
		m->probes.clear(); // mutants dont have probes
		probes.push_back(m);
	}
}


template <class Model>
void MySpecies<Model>::evolveTraits(double dt){
	std::vector<double> traits_res = get_traits();
	std::vector<double> dx(traits_res.size());

	for (int i=0; i<dx.size(); ++i) dx[i] = fitness_gradient[i]*trait_variance[i]*dt;

	for (int i=0; i<dx.size(); ++i) traits_res[i] += dx[i];
	set_traits(traits_res);

	for (auto m : probes){
		std::vector<double> traits_mut = m->get_traits();
		for (int j=0; j<traits_mut.size(); ++j) traits_mut[j] += dx[j];
		m->set_traits(traits_mut);
	}
}

template <class Model>
void MySpecies<Model>::calcFitnessGradient(){
	fitness_gradient.clear();
	for (auto m : probes){
		double grad = log(m->r0/r0) / fg_dx;
		fitness_gradient.push_back(grad);
		//cout << "   " << m->invasion_fitness << " " << spp->invasion_fitness << " " << *spp->fitness_gradient.rbegin() << "\n";
	}
}

template <class Model>
void MySpecies<Model>::print_extra(){
	std::cout << "Resident: " << ((isResident)? "Yes" : "No") << "\n";
	std::cout << "Probes: ";
	for (auto x : probes) std::cout << x << " ";
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