#include "state_restore.h"
using namespace std;

void saveState(Solver *S, string configfilename){
	io::Initializer I(configfilename);
	I.readFile();
	string outdir   = I.get<string>("outDir") + "/" + I.get<string>("exptName");
	string filename = outdir + "/" + I.get<string>("savedState");

	ofstream fout(filename.c_str());
	if (!fout) throw runtime_error("Could not open file for saving state: "+filename);

	fout << "Plant-FATE::state::v1" << '\n';

	fout.close();
}

void restoreState(Solver * S, string configfilename){
	io::Initializer I(configfilename);
	I.readFile();
	string outdir   = I.get<string>("outDir") + "/" + I.get<string>("exptName");
	string filename = outdir + "/" + I.get<string>("savedState");

	ifstream fin(filename.c_str());
	if (!fin) throw runtime_error("Could not open file for restoring state: "+filename);

	string s; fin >> s;  // discard version number

	fin.close();
}

