#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
using namespace std;

#include <solver.h>
#include "pspm_interface.h"
#include "trait_reader.h"

//*********************************************************************************
//**************** NO CHANGES MADE TO ORIGINAL CODE - DO NOT TOUCH ****************
//*********************************************************************************

std::vector <double> myseq(double from, double to, int len){
    std::vector<double> x(len);
    for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
    return x;
}

inline double runif(double rmin=0, double rmax=1){
    double r = double(rand())/RAND_MAX;
    return rmin + (rmax-rmin)*r;
}

vector<double> generateDefaultCohortSchedule(double max_time){

    vector<double> tvec;

    const double multiplier=0.2, min_step_size=1e-5, max_step_size=2;
    
    assert(min_step_size > 0 && "The minimum step size must be greater than zero");
    
    double dt = 0.0, time = 0.0;
    tvec.push_back(time);
    while (time <= max_time) {
        dt = exp2(floor(log2(time * multiplier)));
        time += min(max(dt, min_step_size), max_step_size);
        tvec.push_back(time);
    }

    // Drop the last time; that's not going to be needed:
    if (tvec.size() >=1)     // JAI: added to avoid overflow warning
        tvec.resize(tvec.size() - 1);

    return tvec;
}

class SolverIO{
    public:
    int nspecies;
    Solver * S;

    vector <vector<ofstream>> streams;

    void openStreams(vector<string> varnames, string dir = "."){
        varnames.insert(varnames.begin(), "u");
        varnames.insert(varnames.begin(), "X");
        
        for (int s=0; s < S->species_vec.size(); ++s){
            auto spp = S->species_vec[s];
            vector<ofstream> spp_streams;
            
            for (int i=0; i<varnames.size(); ++i){
                stringstream sout;
                sout << dir << "/species_" << s << "_" << varnames[i] << ".txt";
                cout << sout.str() << endl;
                ofstream fout(sout.str().c_str());
                assert(fout);
                spp_streams.push_back(std::move(fout));
            }
            streams.push_back(std::move(spp_streams));
        }
    }

    void closeStreams(){
        for (int s=0; s<streams.size(); ++s){
            for (int j=0; j<streams[s].size(); ++j){
                streams[s][j].close();
            }
        }
    }

    void writeState(){
        for (int s=0; s < S->species_vec.size(); ++s){
            auto spp = (Species<PSPM_Plant>*)S->species_vec[s];

            for (int i=0; i<streams[s].size(); ++i) streams[s][i] << S->current_time << "\t";

            for (int j=0; j<spp->xsize(); ++j){
                auto& C = spp->getCohort(j);
                int is = 0;
                streams[s][is++] << C.x << "\t";
                streams[s][is++] << C.u << "\t";
                streams[s][is++] << C.geometry.height << "\t";
                streams[s][is++] << C.geometry.lai << "\t";
                streams[s][is++] << C.rates.dmort_dt << "\t";
                streams[s][is++] << C.state.seed_pool << "\t";
                streams[s][is++] << C.rates.rgr << "\t";
                streams[s][is++] << C.res.gpp/C.geometry.crown_area << "\t";

            }
            
            for (int i=0; i<streams[s].size(); ++i) streams[s][i] << endl; //"\n";
        }
    }
};

class CWM{
    public:
    double n_ind;
    double biomass;
    double ba;
    double canopy_area;
    double height;
    double lma;
    double p50;
    double hmat;
    double wd;
    double gs;
    
    vector<double> n_ind_vec;
    vector<double> biomass_vec;
    vector<double> ba_vec;
    vector<double> canopy_area_vec;
    vector<double> height_vec;

    vector<double> lma_vec;
    vector<double> p50_vec;
    vector<double> hmat_vec;
    vector<double> wd_vec;
    
    void update(double t, Solver &S){
        n_ind_vec.clear();
        n_ind_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k)
            n_ind_vec[k] = S.integrate_x([&S,k](int i, double t){
                                              return 1;
                                        }, t, k);
        n_ind = std::accumulate(n_ind_vec.begin(), n_ind_vec.end(), 0.0);

        biomass_vec.clear();
        biomass_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k)
            biomass_vec[k] = S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.get_biomass();
                                        }, t, k);
        biomass = std::accumulate(biomass_vec.begin(), biomass_vec.end(), 0.0);

        ba_vec.clear();
        ba_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k)
            ba_vec[k] = S.integrate_wudx_above([&S,k](int i, double t){
                                              double D = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry.diameter;
                                              return M_PI*D*D/4;
                                        }, t, 0.1, k);
        ba = std::accumulate(ba_vec.begin(), ba_vec.end(), 0.0);

        canopy_area_vec.clear();
        canopy_area_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k)
            canopy_area_vec[k] = S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry;
                                              return p.crown_area;
                                        }, t, k);
        canopy_area = std::accumulate(canopy_area_vec.begin(), canopy_area_vec.end(), 0.0);
        

        height_vec.clear();
        height_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k)
            height_vec[k] = S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.geometry.height;
                                        }, t, k);
                                        
        for (int k=0; k<S.n_species(); ++k) height_vec[k] /= n_ind_vec[k];


        hmat = 0;
        for (int k=0; k<S.n_species(); ++k)
            hmat += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.traits.hmat;
                                        }, t, k);
        hmat /= n_ind;
        hmat_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k) hmat_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.hmat;


        lma = 0;
        for (int k=0; k<S.n_species(); ++k)
            lma += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.traits.lma;
                                        }, t, k);
        lma /= n_ind;
        lma_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k) lma_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.lma;

        wd = 0;
        for (int k=0; k<S.n_species(); ++k)
            wd += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.traits.wood_density;
                                        }, t, k);
        wd /= n_ind;
        wd_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k) wd_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.wood_density;

        p50 = 0;
        for (int k=0; k<S.n_species(); ++k)
            p50 += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.traits.p50_xylem;
                                        }, t, k);
        p50 /= n_ind;
        p50_vec.resize(S.n_species());
        for (int k=0; k<S.n_species(); ++k) p50_vec[k] = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(-1).traits.p50_xylem;

        gs = 0;
        for (int k=0; k<S.n_species(); ++k)
            gs += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.res.gs_avg;
                                        }, t, k);
        gs /= n_ind;
    
    }
};

class EmergentProps{
    public:
    double gpp;
    double npp;
    double resp_auto;
    double trans;
    double lai;
    double leaf_mass;
    double stem_mass;
    double croot_mass;
    double froot_mass;

    void update(double t, Solver &S){
        // gpp
        gpp = 0;
        for (int k=0; k<S.n_species(); ++k)
            gpp += S.integrate_x([&S,k](int i, double t){
                                      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                      return p.res.gpp;
                                  }, t, k);

        // npp
        npp = 0;
        for (int k=0; k<S.n_species(); ++k)
            npp += S.integrate_x([&S,k](int i, double t){
                                      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                      return p.res.npp;
                                  }, t, k);

        // transpiration
        trans = 0;
        for (int k=0; k<S.n_species(); ++k)
            trans += S.integrate_x([&S,k](int i, double t){
                                      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                      return p.res.trans;
                                  }, t, k);


        // autotropic respiration
        resp_auto = 0;
        for (int k=0; k<S.n_species(); ++k)
            resp_auto += S.integrate_x([&S,k](int i, double t){
                                      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                      return p.res.rleaf + p.res.rroot + p.res.rstem;
                                  }, t, k);
        
        // LAI
        lai = 0;
        for (int k=0; k<S.n_species(); ++k)
            lai += S.integrate_x([&S,k](int i, double t){
                                      auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i).geometry;
                                      return p.crown_area*p.lai;
                                }, t, k);

        // Leaf mass
        leaf_mass = 0;
        for (int k=0; k<S.n_species(); ++k)
            leaf_mass += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.geometry.leaf_mass(p.traits);
                                        }, t, k);

        // Wood mass
        stem_mass = 0;
        for (int k=0; k<S.n_species(); ++k)
            stem_mass += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.geometry.stem_mass(p.traits);
                                        }, t, k);
        
        // coarse root mass
        croot_mass = 0;
        for (int k=0; k<S.n_species(); ++k)
            croot_mass += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.geometry.coarse_root_mass(p.traits);
                                        }, t, k);

        // fine root mass
        froot_mass = 0;
        for (int k=0; k<S.n_species(); ++k)
            froot_mass += S.integrate_x([&S,k](int i, double t){
                                              auto& p = (static_cast<Species<PSPM_Plant>*>(S.species_vec[k]))->getCohort(i);
                                              return p.geometry.root_mass(p.traits);
                                        }, t, k);

    }
};

//*********************************************************************************
//************** NEW CODE FOR PATCH STRUCTURE STARTS HERE - CAN EDIT **************
//*********************************************************************************
