#ifndef PHYDRO_TEMPERATURE_DEPENDENCIES_PHOTOSYNTHESIS_H
#define PHYDRO_TEMPERATURE_DEPENDENCIES_PHOTOSYNTHESIS_H

#include <cmath>
#include <stdexcept>

#include "temperature_dependencies_physical.h"  // for calc_patm()

namespace phydro{

// FIXME: use enum class instead?
enum FtempVcmaxJmaxMethod{FV_kattge07, 
                          FV_kumarathunge19, 
                          FV_leuning02};

enum FtempRdMethod{FR_heskel16, 
                   FR_arrhenius, 
                   FR_q10};

enum FtempBrMethod{FB_atkin15, 
                   FB_kumarathunge19};



//-----------------------------------------------------------------------
// Output:   Factor fv to correct for instantaneous temperature response
//           of Vcmax for:
//
//               Vcmax(temp) = fv * Vcmax(25 deg C) 
//
//-----------------------------------------------------------------------
inline float calc_ftemp_arrhenius(float tk, float dha, float tkref = 298.15 ){

  // # Note that the following forms are equivalent:
  // # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
  // # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
  // # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )

  float kR   = 8.3145;     // Universal gas constant, J/mol/K

  float ftemp = exp( dha * (tk - tkref) / (tkref * kR * tk) );

  return ftemp;
}



//-----------------------------------------------------------------------
// Input:    - float, air temperature, deg C (temp)
//           - float, atmospheric pressure, Pa (patm)
// Output:   float, Pa (mmk)
// Features: Returns the temperature & pressure dependent Michaelis-Menten
//           coefficient, K (Pa).
// Ref:      Bernacchi et al. (2001), Improved temperature response 
//           functions for models of Rubisco-limited photosynthesis, 
//           Plant, Cell and Environment, 24, 253--259.
//-----------------------------------------------------------------------
inline float calc_kmm(float tc, float patm ) {
  float dhac   = 79430;      // (J/mol) Activation energy, Bernacchi et al. (2001)
  float dhao   = 36380;      // (J/mol) Activation energy, Bernacchi et al. (2001)
  float kco    = 2.09476e5;  // (ppm) O2 partial pressure, Standard Atmosphere

  //// k25 parameters are not dependent on atmospheric pressure
  float kc25 = 39.97;   // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  float ko25 = 27480;   // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  float tk = tc + 273.15;

  float kc = kc25 * calc_ftemp_arrhenius( tk, dhac );
  float ko = ko25 * calc_ftemp_arrhenius( tk, dhao );

  float po  = kco * (1e-6) * patm;         // O2 partial pressure
  float kmm = kc * (1.0 + po/ko);

  return kmm;
}




//-----------------------------------------------------------------------
// Input:    float, air temperature, degrees C (tc)
// Output:   float, gamma-star, Pa (gammastar)
// Features: Returns the temperature-dependent photorespiratory 
//           compensation point, Gamma star (Pascals), based on constants 
//           derived from Bernacchi et al. (2001) study.
// Ref:      Bernacchi et al. (2001), Improved temperature response 
//           functions for models of Rubisco-limited photosynthesis, 
//           Plant, Cell and Environment, 24, 253--259.
//-----------------------------------------------------------------------
inline float calc_gammastar(float tc, float patm ) {
  float dha    = 37830;       // (J/mol) Activation energy, Bernacchi et al. (2001)
  float gs25_0 = 4.332;       // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  float gammastar25 = gs25_0 * patm / calc_patm(0.0);

  float tk = tc + 273.15;
  float gammastar = gammastar25 * calc_ftemp_arrhenius( tk, dha );

  return gammastar;
}


inline double calc_ftemp_kphio(double tc, bool c4 = false) {
    double ftemp;
    
    if (c4) {
        ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc*tc;
    } 
    else {
        ftemp = 0.352 + 0.022 * tc - 3.4e-4 * tc*tc;
    }
    
    // Avoid negative values
    if (ftemp < 0.0) {
        ftemp = 0.0;
    }
    
    return ftemp;
}


//-----------------------------------------------------------------------
// arguments
// tcleaf: temperature (degrees C)
// tref: is 'to' in Nick's set it to 25 C (=298.15 K in other cals)
//
// function return variable
// fv: temperature response factor, relative to 25 deg C.
//
// Output:   Factor fv to correct for instantaneous temperature response
//           of Vcmax or Jmax for:
//
//               Vcmax(temp) = fv * Vcmax(25 deg C) 
//                Jmax(temp) = fv *  Jmax(25 deg C) 
//
// Ref:      Pascal Schneider et al. (in prep.) Optimal temperature paper
//-----------------------------------------------------------------------
inline double calc_ftemp_inst_vcmax(double tcleaf, double tcgrowth = 1e20, double tcref = 25.0, FtempVcmaxJmaxMethod method_ftemp = FV_kumarathunge19) {
    double Rgas = 8.3145; // Universal gas constant (J/mol/K)
    double tkref = tcref + 273.15; // Convert reference temperature to Kelvin
    double tkleaf = tcleaf + 273.15; // Convert leaf temperature to Kelvin
    double fv;

    if (method_ftemp == FV_kattge07 || method_ftemp == FV_kumarathunge19) {
        // Kattge2007 Parametrization
        double Hd = 200000; // Deactivation energy (J/mol)
        double Ha = 71513; // Activation energy (J/mol)
        double a_ent = 668.39; // Offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
        double b_ent = 1.07; // Slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

        if (method_ftemp == FV_kumarathunge19) {
            // Kumarathunge2019 Implementation:
            // local parameters
            a_ent = 645.13; // Offset of entropy vs. temperature relationship (J/mol/K)
            b_ent = 0.38; // Slope of entropy vs. temperature relationship (J/mol/K^2)
            
            // local variables
            Ha = 42600 + (1140 * tcgrowth); // Acclimation for vcmax
        }

        // Calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin!
        double dent = a_ent - (b_ent * tcgrowth);  //  'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's

        double fva = calc_ftemp_arrhenius(tkleaf, Ha, tkref);
        double fvb = (1 + std::exp((tkref * dent - Hd) / (Rgas * tkref))) / (1 + std::exp((tkleaf * dent - Hd) / (Rgas * tkleaf)));
        fv = fva * fvb;
    } 
    else if (method_ftemp == FV_leuning02) {
        // Ref: Leuning, R. (2002). Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25(9), 1205–1210. https://doi.org/10.1046/j.1365-3040.2002.00898.x
        // Table 2:
        double Ha = 73637;
        double Hd = 149252;
        double Sv = 486;

        double term_1 = 1 + std::exp((Sv * tkref - Hd) / (Rgas * tkref));
        double term_3 = 1 + std::exp((Sv * tkleaf - Hd) / (Rgas * tkleaf));
        double term_2 = std::exp((Ha / (Rgas * tkref)) * (1 - tkref / tkleaf)); // Careful: In Eq. (1) in Leuning et al. (1992), there is a bracket missing in this term!

        fv = term_1 * term_2 / term_3;
    } 
    else {
        throw std::invalid_argument("Invalid method_ftemp:" + method_ftemp);
    }

    return fv;
}


inline double calc_ftemp_inst_jmax(double tcleaf, double tcgrowth, double tchome = 1e20, double tcref = 25.0, FtempVcmaxJmaxMethod method_ftemp = FV_kumarathunge19) {
    double Rgas = 8.3145; // Universal gas constant (J/mol/K)
    double tkref = tcref + 273.15; // Convert reference temperature to Kelvin
    double tkleaf = tcleaf + 273.15; // Convert leaf temperature to Kelvin
    double fv;

    if (method_ftemp == FV_kattge07 || method_ftemp == FV_kumarathunge19) {
        double Hd = 200000; // Deactivation energy (J/mol)
        double Ha = 49884; // Activation energy (J/mol)
        double a_ent = 659.70; // Offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
        double b_ent = 0.75; // Slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

        // Calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin!
        double dent = a_ent - b_ent * tcgrowth;  // 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's

        if (method_ftemp == FV_kumarathunge19) {
            // Kumarathunge2019 Implementation:
            // local parameters
            Ha = 40710; // Activation energy (J/mol)
            a_ent = 658.77; // Offset of entropy vs. temperature relationship (J/mol/K)
            b_ent = 0.84; // Slope of entropy vs. temperature relationship (J/mol/K^2)
            double c_ent = 0.52; // 2nd slope of entropy vs. temperature (J/mol/K^2)

            // Entropy calculation, equations given in Celsius, not in Kelvin
            dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - tchome);
        }

        double fva = calc_ftemp_arrhenius(tkleaf, Ha, tkref);
        double fvb = (1 + std::exp((tkref * dent - Hd) / (Rgas * tkref))) / (1 + std::exp((tkleaf * dent - Hd) / (Rgas * tkleaf)));
        fv = fva * fvb;

    } else if (method_ftemp == FV_leuning02) {
        // Ref: Leuning, R. (2002). Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25(9), 1205–1210. https://doi.org/10.1046/j.1365-3040.2002.00898.x
        // Table 2:
        double Ha = 50300;
        double Hd = 152044;
        double Sv = 495;

        double term_1 = 1 + std::exp((Sv * tkref - Hd) / (Rgas * tkref));
        double term_3 = 1 + std::exp((Sv * tkleaf - Hd) / (Rgas * tkleaf));
        double term_2 = std::exp((Ha / (Rgas * tkref)) * (1 - tkref / tkleaf)); // Careful: In Eq. (1) in Leuning et al. (1992), there is a bracket missing in this term!

        fv = term_1 * term_2 / term_3;
    } else {
        throw std::invalid_argument("Invalid method_ftemp:" + method_ftemp);
    }

    return fv;
}


// Calculate Temperature scaling (f) factor for Rd,
// Rd = f * Rd25
inline double calc_ftemp_inst_rd(double tc_leaf, FtempRdMethod method_rd_scale = FR_heskel16){//, double tc_growth = 1e20, double q10 = 2) {
    // Get temperature scaling for Rd:

    double f = 1.0; // Scaling factor for Rd

    if (method_rd_scale == FR_heskel16) {
        // Heskel et al. (2016) temperature scaling
        double apar = 0.1012;
        double bpar = 0.0005;
        f = std::exp(apar * (tc_leaf - 25.0) - bpar * (tc_leaf*tc_leaf - 25.0*25.0));
    }
    else if (method_rd_scale == FR_arrhenius) {
        // Arrhenius temperature scaling
        double dha = 20700; // Activation energy taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
        f = calc_ftemp_arrhenius(tc_leaf + 273.15, dha); // Convert temperature to Kelvin and call calc_ftemp_arrh function
    }
    else if (method_rd_scale == FR_q10) {
        // Q10 temperature scaling according to Tjoelker et al. (2001)
        f = std::pow(3.22 - 0.046 * tc_leaf, (tc_leaf - 25.0)) / 10;
    }
    else {
        throw std::invalid_argument("Invalid method_rd_scale:" + method_rd_scale);
    }

    return f;
}

// Ratio of Rd to Vcmax at 25 degC
// Rd25 = brd25 * Vcmax25 
inline double calc_brd25(FtempBrMethod method_rd25 = FB_atkin15, double tc_growth = 25.0) {
    double rd_to_vcmax;

    if (method_rd25 == FB_atkin15) {
        rd_to_vcmax = 0.015; // Ratio of Rdark to Vcmax25, Atkin et al., 2015 for C3 herbaceous
    }
    else if (method_rd25 == FB_kumarathunge19) {
        rd_to_vcmax = 0.0360 - 0.0010 * tc_growth; // Acclimated rd_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
    }
    else {
        throw std::invalid_argument("Invalid method_rd25:" + method_rd25);
    }

    return rd_to_vcmax;
}



//-----------------------------------------------------------------------
// arguments
// tcleaf: temperature (degrees C)
// tref: is 'to' in Nick's set it to 25 C (=298.15 K in other cals)
//
// function return variable
// fv: temperature response factor, relative to 25 deg C.
//
// Output:   Factor fv to correct for instantaneous temperature response
//           of Vcmax for:
//
//               Vcmax(temp) = fv * Vcmax(25 deg C) 
//
// Ref:      Wang Han et al. (in prep.)
//-----------------------------------------------------------------------
inline float calc_ftemp_inst_vcmax_WangEtAl(float tcleaf, float tcgrowth, float tcref = 25.0 ){
  // loal parameters
  float Ha    = 71513;  // activation energy (J/mol)
  float Hd    = 200000; // deactivation energy (J/mol)
  float Rgas  = 8.3145; // universal gas constant (J/mol/K)
  float a_ent = 668.39; // offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
  float b_ent = 1.07;   // slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
  
  float tkref = tcref + 273.15;  // to Kelvin

  // conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C. 
  float tkleaf = tcleaf + 273.15;

  // calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
  float dent = a_ent - b_ent * tcgrowth;   // 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
  float fva = calc_ftemp_arrhenius( tkleaf, Ha, tkref );
  float fvb = (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) );
  float fv  = fva * fvb;

  return fv;
}

inline float calc_ftemp_vcmax_bernacchi(double tc){
  double dha = 65.33e3;
  double c = 26.35;
  double R = 8.31;
  return exp(c - dha/R/(tc+273.16));
}

//-----------------------------------------------------------------------
// arguments:
// tc                  // temperature (degrees C)
// function return variable:
// fr                  // temperature response factor, relative to 25 deg C.
// Output:   Factor fr to correct for instantaneous temperature response
//           of Rd (dark respiration) for:
//
//               Rd(temp) = fr * Rd(25 deg C) 
//
// Ref:      Heskel et al. (2016) used by Wang Han et al. (in prep.)
//-----------------------------------------------------------------------
inline float calc_ftemp_inst_rd_heskel_only(float tc){
  // loal parameters
  float apar = 0.1012;
  float bpar = 0.0005;
  float tk25  = 298.15; // 25 deg C in Kelvin

  // conversion of temperature to Kelvin
  float tk = tc + 273.15;

  float fr = exp( apar * (tc - 25.0) - bpar * (tc*tc - 25.0*25.0) );
 
  return fr; 
}


} // phydro


#endif