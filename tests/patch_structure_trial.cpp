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

class Patch{
    public:
    
    
    void initPatch(PSPM_Dynamic_Environment &E, TraitsReader &Tr, int nspp, double y0, int patchNo){
        Solver S(SOLVER_IFMU, "rk45ck");
        S.control.ode_ifmu_stepsize = 0.0833333;
        S.control.ifmu_centered_grids = false; //true;
        S.use_log_densities = true;
        S.setEnvironment(&E);
        
        for (int i=0; i<nspp; ++i){
            PSPM_Plant p1;
            p1.initParamsFromFile("tests/params/p.ini");
            p1.traits.species_name = Tr.species[i].species_name;
            p1.traits.lma = Tr.species[i].lma;
            p1.traits.wood_density = Tr.species[i].wood_density;
            p1.traits.hmat = Tr.species[i].hmat;
            p1.traits.p50_xylem = Tr.species[i].p50_xylem; // runif(-3.5,-0.5);
            
            p1.coordinateTraits();

            ((plant::Plant*)&p1)->print();
            
            //p1.geometry.set_lai(p1.par.lai0); // these are automatically set by init_state() in pspm_interface
            p1.set_size(0.01);
            
            Species<PSPM_Plant>* spp = new Species<PSPM_Plant>(p1);

            S.addSpecies(30, 0.01, 10, true, spp, 3, 1e-3);
            
            //    S.addSpecies(vector<double>(1, p1.geometry.get_size()), &spp, 3, 1);
            //S.get_species(0)->set_bfin_is_u0in(true);    // say that input_birth_flux is u0
        }
        S.resetState(y0);
        S.initialize();
        
        for (auto spp : S.species_vec) spp->setU(0, 1);
        S.copyCohortsToState();

        S.print();
    //  S.control.update_cohorts = false;

        SolverIO sio;
        sio.S = &S;
        sio.openStreams({"height", "lai", "mort", "seeds", "g", "gpp"}, string(out_dir+"/Patch"+patchNo);


        //-------------------------------------//
        // CREATING FILES AND ADDING DATA
        //-------------------------------------//


        string patchpath = string(out_dir + "/Patch"+ patchNo);
        string command = "mkdir -p " + patchpath;
        int sysresult;
        sysresult = system(command.c_str());
        //ADD PATCH NUMBER FOR EACH FILE
        //  ofstream fout("fmu_PlantFATE.txt");
        ofstream fzst(string(patchpath + "/z_star.txt").c_str());
        assert(fzst);
        ofstream fco(string(patchpath + "/canopy_openness.txt").c_str());
        ofstream fseed(string(patchpath + "/seeds.txt").c_str());
        ofstream fabase(string(patchpath + "/basal_area.txt").c_str());
        //	ofstream flai(string(out_dir + "/LAI.txt").c_str());
        //	ofstream fcwmt(string(out_dir + "/cwmt.txt").c_str());
        ofstream foutd(string(patchpath + "/" + I.get<string>("emgProps")).c_str());
        ofstream fouty(string(patchpath + "/" + I.get<string>("cwmAvg")).c_str());
        ofstream fouty_spp(string(patchpath + "/" + I.get<string>("cwmperSpecies")).c_str());


        foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\n";
        fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
        fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
        double t_clear = 500;
        // t is years since 2000-01-01
        double y0, yf;
        y0 = I.getScalar("year0");
        yf = I.getScalar("yearf");

        for (double t=y0; t <= yf; t=t+1) {
            cout << "t = " << t << endl; //"\t";
            S.step_to(t, after_step);
            
            CWM cwm;
            EmergentProps props; 
            
            cwm.update(t, S);
            props.update(t, S);
            
            foutd << int(t) << "\t"
                << (t-int(t))*365 << "\t"
                << props.gpp*0.5/365*1000 << "\t"
                << props.npp*0.5/365*1000 << "\t"
                << props.resp_auto*0.5/365*1000 << "\t"  // gC/m2/d
                << props.leaf_mass*1000*0.5 << "\t"     
                << props.stem_mass*1000*0.5 << "\t"
                << props.croot_mass*1000*0.5 << "\t"
                << props.froot_mass*1000*0.5 << "\t"
                << (props.croot_mass+props.froot_mass)*1000*0.5 << "\t" // gC/m2
                << cwm.gs << "\t"
                << props.trans/365 << "\t"   // kg/m2/yr --> 1e-3 m3/m2/yr --> 1e-3*1e3 mm/yr --> 1/365 mm/day  
                << props.lai << endl;
            
            fouty << int(t) << "\t"
                << -9999  << "\t"
                << cwm.n_ind << "\t"
                << -9999  << "\t"
                << cwm.height  << "\t"
                << cwm.hmat  << "\t"
                << cwm.canopy_area  << "\t"   // m2/m2
                << cwm.ba  << "\t"            // m2/m2
                << cwm.biomass  << "\t"       // kg/m2
                << cwm.wd  << "\t"
                << -9999  << "\t"
                << 1/cwm.lma  << "\t"
                << cwm.p50  << endl;
            
            for (int k=0; k<S.species_vec.size(); ++k){
                fouty_spp 
                    << int(t) << "\t"
                    << k  << "\t"
                    << cwm.n_ind_vec[k] << "\t"
                    << -9999  << "\t"
                    << cwm.height_vec[k]  << "\t"
                    << cwm.hmat_vec[k]  << "\t"
                    << cwm.canopy_area_vec[k]  << "\t"   // m2/m2
                    << cwm.ba_vec[k]  << "\t"            // m2/m2
                    << cwm.biomass_vec[k]  << "\t"       // kg/m2
                    << cwm.wd_vec[k]  << "\t"
                    << -9999  << "\t"
                    << 1/cwm.lma_vec[k]  << "\t"
                    << cwm.p50_vec[k]  << "\n";
            }

            fseed << t << "\t";
            for (int i=0; i<S.n_species(); ++i) fseed << seeds_hist[i].get() << "\t";
            fseed << "\n";

            fabase << t << "\t";
            for (int i=0; i<S.n_species(); ++i) fabase << cwm.ba_vec[i] << "\t";
            fabase << "\n";

            		fzst << t << "\t";
            for (auto z : E.z_star) fzst << z << "\t";
            fzst << endl;
            
            fco << t << "\t";
            for (auto z : E.canopy_openness) fco << z << "\t";
            fco << endl;
                
            sio.writeState();
        
            // clear patch after 50 year	
            if (t >= t_clear){
                for (auto spp : S.species_vec){
                    for (int i=0; i<spp->xsize(); ++i){
                        auto& p = (static_cast<Species<PSPM_Plant>*>(spp))->getCohort(i);
                        p.geometry.lai = p.par.lai0;
                        double u_new = 0; //spp->getU(i) * double(rand())/RAND_MAX;
                        spp->setU(i, u_new);
                    }
                }
                S.copyCohortsToState();
                double t_int = -log(double(rand())/RAND_MAX) * I.getScalar("T_return");;
                t_clear = t + fmin(t_int, 1000);
            }
            
            fco.flush();
            fseed.flush();
            fzst.flush();
            fabase.flush();
    //		flai.flush();
    //		fcwmt.flush();

        }
        

        fco.close();
        fseed.close();
        fzst.close();
        fabase.close();
    //	flai.close();
    //	fcwmt.close();
        foutd.close();
        fouty.close();
        // free memory associated
        for (auto s : S.species_vec) delete static_cast<Species<PSPM_Plant>*>(s); 
           
    }
    
    void seed_thing(){
        
    }
    
//     void create_files(int patchNo){
//         string patchpath = string(out_dir + "/Patch"+ patchNo);
//         string command = "mkdir -p " + patchpath;
//         int sysresult;
//         sysresult = system(command.c_str());
//         //ADD PATCH NUMBER FOR EACH FILE
//         //  ofstream fout("fmu_PlantFATE.txt");
//         ofstream fzst(string(patchpath + "/z_star.txt").c_str());
//         assert(fzst);
//         ofstream fco(string(patchpath + "/canopy_openness.txt").c_str());
//         ofstream fseed(string(patchpath + "/seeds.txt").c_str());
//         ofstream fabase(string(patchpath + "/basal_area.txt").c_str());
//         //	ofstream flai(string(out_dir + "/LAI.txt").c_str());
//         //	ofstream fcwmt(string(out_dir + "/cwmt.txt").c_str());
//         ofstream foutd(string(patchpath + "/" + I.get<string>("emgProps")).c_str());
//         ofstream fouty(string(patchpath + "/" + I.get<string>("cwmAvg")).c_str());
//         ofstream fouty_spp(string(patchpath + "/" + I.get<string>("cwmperSpecies")).c_str());


//         foutd << "YEAR\tDOY\tGPP\tNPP\tRAU\tCL\tCW\tCCR\tCFR\tCR\tGS\tET\tLAI\n";
//         fouty << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
//         fouty_spp << "YEAR\tPID\tDE\tOC\tPH\tMH\tCA\tBA\tTB\tWD\tMO\tSLA\tP50\n";
//         double t_clear = 500;
//         // t is years since 2000-01-01
//         double y0, yf;
//         y0 = I.getScalar("year0");
//         yf = I.getScalar("yearf");

//         for (double t=y0; t <= yf; t=t+1) {
//             cout << "t = " << t << endl; //"\t";
//             S.step_to(t, after_step);
            
//             CWM cwm;
//             EmergentProps props; 
            
//             cwm.update(t, S);
//             props.update(t, S);
            
//             foutd << int(t) << "\t"
//                 << (t-int(t))*365 << "\t"
//                 << props.gpp*0.5/365*1000 << "\t"
//                 << props.npp*0.5/365*1000 << "\t"
//                 << props.resp_auto*0.5/365*1000 << "\t"  // gC/m2/d
//                 << props.leaf_mass*1000*0.5 << "\t"     
//                 << props.stem_mass*1000*0.5 << "\t"
//                 << props.croot_mass*1000*0.5 << "\t"
//                 << props.froot_mass*1000*0.5 << "\t"
//                 << (props.croot_mass+props.froot_mass)*1000*0.5 << "\t" // gC/m2
//                 << cwm.gs << "\t"
//                 << props.trans/365 << "\t"   // kg/m2/yr --> 1e-3 m3/m2/yr --> 1e-3*1e3 mm/yr --> 1/365 mm/day  
//                 << props.lai << endl;
            
//             fouty << int(t) << "\t"
//                 << -9999  << "\t"
//                 << cwm.n_ind << "\t"
//                 << -9999  << "\t"
//                 << cwm.height  << "\t"
//                 << cwm.hmat  << "\t"
//                 << cwm.canopy_area  << "\t"   // m2/m2
//                 << cwm.ba  << "\t"            // m2/m2
//                 << cwm.biomass  << "\t"       // kg/m2
//                 << cwm.wd  << "\t"
//                 << -9999  << "\t"
//                 << 1/cwm.lma  << "\t"
//                 << cwm.p50  << endl;
            
//             for (int k=0; k<S.species_vec.size(); ++k){
//                 fouty_spp 
//                     << int(t) << "\t"
//                     << k  << "\t"
//                     << cwm.n_ind_vec[k] << "\t"
//                     << -9999  << "\t"
//                     << cwm.height_vec[k]  << "\t"
//                     << cwm.hmat_vec[k]  << "\t"
//                     << cwm.canopy_area_vec[k]  << "\t"   // m2/m2
//                     << cwm.ba_vec[k]  << "\t"            // m2/m2
//                     << cwm.biomass_vec[k]  << "\t"       // kg/m2
//                     << cwm.wd_vec[k]  << "\t"
//                     << -9999  << "\t"
//                     << 1/cwm.lma_vec[k]  << "\t"
//                     << cwm.p50_vec[k]  << "\n";
//             }

//             fseed << t << "\t";
//             for (int i=0; i<S.n_species(); ++i) fseed << seeds_hist[i].get() << "\t";
//             fseed << "\n";

//             fabase << t << "\t";
//             for (int i=0; i<S.n_species(); ++i) fabase << cwm.ba_vec[i] << "\t";
//             fabase << "\n";

//             		fzst << t << "\t";
//             for (auto z : E.z_star) fzst << z << "\t";
//             fzst << endl;
            
//             fco << t << "\t";
//             for (auto z : E.canopy_openness) fco << z << "\t";
//             fco << endl;
                
//             sio.writeState();
        
//             // clear patch after 50 year	
//             if (t >= t_clear){
//                 for (auto spp : S.species_vec){
//                     for (int i=0; i<spp->xsize(); ++i){
//                         auto& p = (static_cast<Species<PSPM_Plant>*>(spp))->getCohort(i);
//                         p.geometry.lai = p.par.lai0;
//                         double u_new = 0; //spp->getU(i) * double(rand())/RAND_MAX;
//                         spp->setU(i, u_new);
//                     }
//                 }
//                 S.copyCohortsToState();
//                 double t_int = -log(double(rand())/RAND_MAX) * I.getScalar("T_return");;
//                 t_clear = t + fmin(t_int, 1000);
//             }
            
//             fco.flush();
//             fseed.flush();
//             fzst.flush();
//             fabase.flush();
//     //		flai.flush();
//     //		fcwmt.flush();

//         }
        

//         fco.close();
//         fseed.close();
//         fzst.close();
//         fabase.close();
//     //	flai.close();
//     //	fcwmt.close();
//         foutd.close();
//         fouty.close();
//         // free memory associated
//         for (auto s : S.species_vec) delete static_cast<Species<PSPM_Plant>*>(s); 

// };

int main(){
    string paramsFile = "tests/params/p.ini";
    io::Initializer I(paramsFile);
    I.readFile();
    string out_dir = I.get<string>("outDir") + "/" + I.get<string>("exptName");
    string command = "mkdir -p " + out_dir;
    string command2 = "cp tests/params/p.ini " + out_dir;
    int sysresult;
    sysresult = system(command.c_str());
    sysresult = system(command2.c_str());
    
    int nspp = I.getScalar("nSpecies");
    double y0 = I.getScalar("year0");
    
    PSPM_Dynamic_Environment E;
    E.metFile = I.get<string>("metFile");
    E.co2File = I.get<string>("co2File");
    E.init();
    E.print(0);
    E.use_ppa = true;
    E.update_met = true;
    E.update_co2 = true;
    
    Solver S(SOLVER_IFMU, "rk45ck");
    S.control.ode_ifmu_stepsize = 0.0833333;
    S.control.ifmu_centered_grids = false; //true;
    S.use_log_densities = true;
    S.setEnvironment(&E);
    
    TraitsReader Tr;
    Tr.readFromFile(I.get<string>("traitsFile"));
    Tr.print();
    
    int npatch = I.getScalar("nPatches");
    for(int i=0; i<npatch; ++i){
        Patch pa1;
        // pa1.create_files(i);
        pa1.initPatch(&E, &Tr, nspp, y0, i);
    }
    
    double T_seed_rain_avg = I.getScalar("T_seed_rain_avg");
    
    //seed things
    
    //open file for each patch
    
    
    // ofstream fco[npatch];
    
    // fco[i](string(out_dir + "/canopy_openness.txt").c_str());
    
}
// END OF FILE
