#include "state_restore.h"
#include <filesystem>
using namespace std;

void saveState(Solver *S, string state_outfile, string config_outfile, string params_file){

	cout << "Saving state to: " << state_outfile << '\n';
	cout << "Saving config to: " << config_outfile << '\n';

	// save config file and get filename of saved state
	// string command = "cp " + params_file + " " + config_outfile;
	// int sysresult = system(command.c_str());
	if (std::filesystem::exists(config_outfile)) std::filesystem::remove(config_outfile); // use this because the overwrite flag in below command does not work!
	std::filesystem::copy(params_file, config_outfile, std::filesystem::copy_options::overwrite_existing);

	// open file for writing state
	ofstream fout(state_outfile.c_str());
	if (!fout) throw runtime_error("Could not open file for saving state: "+state_outfile);

	fout << setprecision(12);
	// core state writing
	fout << "Plant-FATE::state::v1" << '\n';

	// write species names vector
	fout << S->species_vec.size() << " | ";
	for (auto s : S->species_vec){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(s);
		fout << std::quoted(spp->species_name) << ' ';
	}
	fout << '\n';

	// write species associations (list of probes)
	for (auto s : S->species_vec){
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(s);
		fout << std::quoted(spp->species_name) << ' ';
		fout << spp->probes.size() << " | "; 
		for (auto p : spp->probes) fout << std::quoted(p->species_name) << ' ';
		fout << '\n';
	}

	// save Solver
	S->save(fout);

	fout.close();
}


void restoreState(Solver * S, string state_infile, string config_infile){

	cout << "Restoring state from: " << state_infile << '\n';
	cout << "Restoring config from: " << config_infile << '\n';

	ifstream fin(state_infile.c_str());
	if (!fin) throw runtime_error("Could not open file for restoring state: "+state_infile);

	string s; fin >> s;  // discard version number

	// Read species associations (probes)
	vector<string> spp_names;   
	map<string, int> indices; // name --> index --- just a map that tells at which index in the species vector the species is located
	vector<vector<int>> probe_lists;
	int n; fin >> n >> s; // number of species to read, skip s = " | "
	spp_names.resize(n); 
	probe_lists.resize(n);
	for (int i=0; i<n; ++i) fin >> std::quoted(spp_names[i]);

	// map indices
	for (int i=0; i<n; ++i){
		indices[spp_names[i]] = i;
	}

	// Species read: 
	cout << "Species read:\n";
	for (int i=0; i<spp_names.size(); ++i) cout << i << " " << indices[spp_names[i]] << " " << spp_names[i] << "\n";

	for (string s : spp_names){
		string r_name; vector<string> probes_list;
		// read resident name
		fin >> std::quoted(r_name);
		assert(r_name == s);
		// read probe names
		fin >> n >> s; // s has " | " 
		cout << spp_names[indices[r_name]] << " --> " << n << " " << s;
		for (int i=0; i<n; ++i){
			fin >> std::quoted(s);
			cout << spp_names[indices[s]] << " ";
			probe_lists[indices[r_name]].push_back(indices[s]);
		}
		cout << '\n';
	}

	PSPM_Plant p;
	vector<Species_Base*> spp_proto;
	for (int i=0; i<spp_names.size(); ++i){
		auto spp = new AdaptiveSpecies<PSPM_Plant>(p);
		spp->configfile_for_restore = config_infile; // set config file that AdaptiveSpecies will use to restore cohorts FIXME: Maybe this can be prevented by simply writing parameters only once? - but configfile will be needed even after continuation for inserting new species (immigration)
		spp_proto.push_back(static_cast<Species_Base*>(spp));
	}

	for (auto spp : spp_proto){
		static_cast<Species_Base*>(spp)->print();
	}

	// restore solver
	S->restore(fin, spp_proto);

	// recreate species associations
	for (string s : spp_names){
		int res_id = indices[s];
		auto spp = static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[res_id]);
		for (int probe_id : probe_lists[res_id]){
			spp->probes.push_back(static_cast<AdaptiveSpecies<PSPM_Plant>*>(S->species_vec[probe_id]));
		}	
	}

	fin.close();
}

