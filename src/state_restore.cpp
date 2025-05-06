#include "state_restore.h"
#include "plantfate_patch.h"
#include <filesystem>
using namespace std;

namespace pfate{

void saveState(Patch& P, string state_outfile, string config_outfile, string params_file){

	cout << "Saving state to: " << state_outfile << '\n';
	cout << "Saving config to: " << config_outfile << '\n';

	// save config file and get filename of saved state
	// string command = "cp " + params_file + " " + config_outfile;
	// int sysresult = system(command.c_str());
	if (std::filesystem::exists(config_outfile)) std::filesystem::remove(config_outfile); // use this because the overwrite flag in below command does not work!
	std::filesystem::copy(params_file, config_outfile, std::filesystem::copy_options::overwrite_existing);

	// open file for writing state
	ofstream fout(state_outfile.c_str());
	if (!fout) throw runtime_error("Could not open file for saving state: " + state_outfile);

	fout << setprecision(12);

	// core state writing
	fout << "Plant-FATE::state::v1" << '\n';

	P.save(fout);

	fout.close();
}


void restoreState(Patch& P, string state_infile, string config_infile){

	cout << "Restoring state from: " << state_infile << '\n';
	cout << "Restoring config from: " << config_infile << '\n';

	ifstream fin(state_infile.c_str());
	if (!fin) {
				cerr << "Error opening file: " << state_infile
             << std::endl;

        // Check for specific error conditions
        if (fin.bad()) {
            cerr << "Fatal error: badbit is set." << endl;
        }

        if (fin.fail()) {
            // Print a more detailed error message using
            // strerror
            cerr << "Error details: " << strerror(errno)
                 << endl;
        }
		throw runtime_error("Could not open file for restoring state: " + state_infile);
	}
	string s; fin >> s;  // discard version number

	P.restore(fin);

	fin.close();
}

} // namespace pfate