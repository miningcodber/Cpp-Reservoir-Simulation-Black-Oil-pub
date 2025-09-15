#include <FluidProperties.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

//---- Measurment Units ----//
/*
    Pressure            ---->> Pascal (Pa)
    Temperature         ---->> Kelvin (K)
    Molecular Weight    ---->>        (gr/mol)
    Specific Gravity    ---->> Dimensionless
    Compressibility     ---->>        (1/Pa)
    Gas Constant        ---->> 8.314  (J/mol*K)
    Viscosity           ---->>        (Cp)
    Density             ---->>        (g/cm^3)
*/

//Namespace for Unit Conversions
namespace UnitConversion{
    //-------------------Common Conversions-------------------
    inline double Pa_to_Psia(double presusre_Pa){
        return presusre_Pa / 6894.76;
    }

    inline double Kelvin_to_Rankine(double temperature_K){
        return temperature_K * 1.8;
    }

    inline double Kelvin_to_Fahrenheit(double temperature_K){
        return (temperature_K - 273.15) * 9.0 / 5.0 + 32.0;
    }

    inline double gmol_to_lbmol(double molWeight_gmol){
        return molWeight_gmol / 453.59237;
    }

    inline double g_per_cm3_to_kg_per_m3(double rho){
        return rho * 1000.0;
    }

    inline double API_from_SG(double SG){
        return (141.5 / SG) - 131.5;
    }

    //-------------------1. Gas PVT-------------------
    inline void convertGasPVTInput_SI(double& pressure_Pa, double& temperature_K, double& Pc_Pa, double& Tc_K){
        pressure_Pa = Pa_to_Psia(pressure_Pa);
        Pc_Pa = Pa_to_Psia(Pc_Pa);
        temperature_K = Kelvin_to_Rankine(temperature_K);
        Tc_K = Kelvin_to_Rankine(Tc_K);
    }

    //-------------------2. Gas Viscosity-------------------
    inline void convertGasViscosityInput_SI(double& pressure_Pa, double& temperature_K, double& Pc_Pa, double& Tc_K, double& Mg_gmol){
        convertGasPVTInput_SI(pressure_Pa, temperature_K, Pc_Pa, Tc_K);
        Mg_gmol = Mg_gmol;
    }

    //-------------------3. Oil PVT-------------------
    inline void convertOilPVTInput_SI(double& pressure_Pa, double& temperature_K, double& Psep_Pa, double& Tsep_K, double& API, double Mg_gmol, double SG){
        pressure_Pa = Pa_to_Psia(pressure_Pa);
        temperature_K = Kelvin_to_Fahrenheit(temperature_K);
        Psep_Pa = Pa_to_Psia(Psep_Pa);
        Tsep_K = Kelvin_to_Fahrenheit(Tsep_K);
        API = API_from_SG(SG);
    }

    //-------------------4. Oil PVT-------------------
    inline void convertOilViscosityInput_SI(double& pressure_Pa, double& temperature_K, double& Pb_Pa){
        pressure_Pa = Pa_to_Psia(pressure_Pa);
        Pb_Pa = Pa_to_Psia(Pb_Pa);
        temperature_K = Kelvin_to_Rankine(temperature_K);
    }
}


//Namespace for Correlations
namespace Correlations{
//-------------------1. Gas PVT Correlations ---- BWD - 8constant-------------------

/*
    Pressure       ---->> Psia   
    Temperature    ---->> Rankine (R)
*/
    struct BWRConstants{
        double  A1 = 0.31506237, 
                A2 = -1.04670990, 
                A3 = -0.57832729, 
                A4 = 0.53530771, 
                A5 = -0.61232032, 
                A6 = -0.10488813, 
                A7 = 0.68157001, 
                A8 = 0.68446549;
    };
    //Gas Compressibility
    inline double computeZ_BWR8(double pressure, double temperature, 
                                double Pc, double Tc, 
                                const BWRConstants& A, double tolerance = 1e-6, 
                                double maxIter = 100){
        
        double Pr = pressure / Pc;
        double Tr = temperature / Tc;

        double Z = 1.0; //Initial guess

        //Fixed Poin iteration (Newton Raphson for Z calculation), as Z is implicit Z = f(Z) (appears both on rho_r and Z equations)
        for (int i = 0; i < maxIter; ++i){
            double rho_r = (0.27 * Pr) / (Z * Tr); //Reduced density

            double term1 = (A.A1 + (A.A2 / Tr) + (A.A3 / (Tr * Tr * Tr))) * rho_r;
            double term2 = (A.A4 + (A.A5 / Tr)) * (rho_r * rho_r);
            double term3 = (A.A5 * A.A6) * (std::pow(rho_r, 5) / Tr);
            double term4 = A.A7 * ((rho_r * rho_r) / (Tr * Tr * Tr)) * (1.0 + (A.A8 * rho_r * rho_r)) * std::exp(-A.A8 * rho_r * rho_r);

            double Z_new = 1.0 + term1 + term2 + term3 + term4;
            if(std::abs(Z_new - Z) < tolerance){
                return Z_new;
            }

            Z = Z_new;
        }
        throw std::runtime_error("BWR-8 Z-factor did not converge.");
    }

    //Gas reduced density calculation
    inline double computeRhor_BWR(double Pc, double Tc, double pressure, double temperature, const BWRConstants& A, double tolerance = 1e-6, 
                                  double maxIter = 100){
        double Z = computeZ_BWR8(pressure, temperature, Pc, Tc, A, tolerance, maxIter);

        double Pr = pressure / Pc;
        double Tr = temperature / Tc;
        double rho_r = (0.27 * Pr) / (Z * Tr);

        return rho_r;
    }

    //Gas Fromation Volume Factor
    inline double computeBg_BWR8(double pressure, double temperature, double Pc, double Tc, const BWRConstants& A){
        double tolerance = 1e-6;
        double maxIter = 100;
        double R = 8.3144621; // J/mol * K

        double Z = computeZ_BWR8(pressure, temperature, Pc, Tc, A, tolerance, maxIter);
        double Bg = (Z * R * temperature) / pressure;

        return Bg;
    }

    //Compute isothermal gas compressibility Cr [1/Pa]
    inline double computeCr_BWR8(double pressure, double temperature, 
                                 double Pc, double Tc, 
                                 const BWRConstants& A){
        double Pr = pressure / Pc;
        double Tr = temperature / Tc;
        double tolerance = 1e-6;
        double maxIter = 100;
        double Z = computeZ_BWR8(pressure, temperature, Pc, Tc, A, tolerance, maxIter);
        double rho_r = computeRhor_BWR(Pc, Tc, pressure, temperature, A, tolerance, maxIter);

        //Calculating (θZ/θPr)Tr
        double term1 = (A.A1 + (A.A2/Tr) + (A.A3/(Tr * Tr * Tr)));
        double term2 = 2 * (A.A4 + (A.A5 / Tr)) * rho_r;
        double term3 = 5 * (A.A5 * A.A6) * ((std::pow(rho_r, 4)) / Tr);
        double term4 = 2 * A.A7 * (rho_r / (Tr * Tr * Tr)) * ((1 + (A.A8 * rho_r * rho_r) - (A.A8 * A.A8 * std::pow(rho_r, 4)))) * std::exp(-A.A8 * rho_r * rho_r);
        double dZdPr = term1 + term2 + term3 + term4;
                                    
        //Compressibility Cr = C * Pc
        double Cr = (1/Pr) - (0.27 / (Z * Z * Tr)) * (dZdPr / (1 + (rho_r / Z) * dZdPr));

        double C = Cr / Pc;

        return C;
    }

//-------------------2. Gas Viscosity Correlations -- Lee - Gonzalez - Eakin-------------------

    /*
        Pressure        ---->> (Psia)
        Temperature     ---->> (R)
        Gas Viscosity   ---->> (cp)
        gas density     ---->> (g/cm^3)
    */

    //Gas Viscocity
    inline double computeRho_G(double Pc, double Tc, double pressure, double temperature, double Mg, const BWRConstants& A){
        double tolerance = 1e-6;
        double maxIter = 100;
        double R = 8.3144621; // J/mol * K
        double rhoR = computeRhor_BWR(Pc, Tc, pressure, temperature, A, tolerance, maxIter);
        double rhoG = rhoR * ((Pc * Mg) / (R * Tc));

        return rhoG;
    };

    inline double computeMuG_LGE(double Mg, double pressure, double temperature, double Pc, double Tc, const BWRConstants& A){ //Mg is the Molecular weight of as g/mol
        double A1 = ((9.379 + 0.1607 * Mg) * std::pow(temperature, 1.5)) / (209.2 + (19.26 * Mg) + temperature);
        double A2 = 3.448 + (986.4/temperature) + (0.01009 * Mg);
        double A3 = 2.447 - (0.2224 * A2);

        double tolerance = 1e-6;
        double maxIter = 100;
        double R = 8.3144621; // J/mol * K

        //Calculating rho_g with precalculated rho_r: rho_g = rho_r * ((Pc * Mg)/(R * Tc)) in S.I.

        double Z = computeZ_BWR8(pressure, temperature, Pc, Tc, A, tolerance, maxIter);
        double rho_r = computeRhor_BWR(Pc, Tc, pressure, temperature, A, tolerance, maxIter);
        double rho_g = rho_r * ((Pc * Mg) / (R * Tc));
        double mu_g = A1 * (std::pow(10, -4) * std::exp(A2 * std::pow(rho_g, A3)));

        return mu_g;
    }


//-------------------3. Oil PVT Correlations -- Vasquez & Beggs : two ranges API<30 & API>30-------------------

    /*
        Pressure        ----> (Psia)
        Temperature     ----> (F)
        Compressibility ----> (1/Psia)
        Rsb             ----> (scf/STB)
    */

    struct VBConstants_U30{ //coefficients for γAPI vaalues of under thirty
        double A1 = 4.677 * std::pow(10, -4);
        double A2 = 1.751 * std::pow(10, -5);
        double A3 = -1.811 * std::pow(10, -8);
        double C1 = 0.0362;
        double C2 = 1.0937;
        double C3 = 25.7245;
    };

    struct VBConstants_O30{ //coefficients for γAPI vaalues of over thirty
        double A1 = 4.670 * std::pow(10, -4);
        double A2 = 1.100 * std::pow(10, -5);
        double A3 = 1.377 * std::pow(10, -9);
        double C1 = 0.0178;
        double C2 = 1.1870;
        double C3 = 23.931;
    };

    //Solution Gas Oil Ratio
    inline double Gamma_g(double Mg){
        double gamma_ = Mg / 28.97;

        return gamma_;
    }
    inline double Gamma_g100(double Mg, double API, double Psep, double Tsep){
        double gamma_g = Mg; //removed Gamma_g function implementation
        double gamma_g100 = gamma_g * (1.0 + (5.912 * std::pow(10, -5) * API * Tsep * std::log(Psep / 114.67)));

        return gamma_g100;
    }

    inline double Rso(const VBConstants_U30& U, const VBConstants_O30& O, 
               double pressure, double temperature, double API, double Mg, double Psep, double Tsep, double Pb){
        double rso;
        double gamma_100 = Gamma_g100( Mg,  API,  Psep,  Tsep);
        if(API < 30){
            if(pressure < Pb){
                rso = U.C1 * gamma_100 * std::pow(pressure, U.C2) * (U.C3 * (API / (temperature + 459.67)));
            }else{
                rso = U.C1 * gamma_100 * std::pow(Pb, U.C2) * (U.C3 * (API / (temperature + 459.67)));
            }
        }else{
            if(pressure < Pb){
                rso = O.C1 * gamma_100 * std::pow(pressure, O.C2) * (O.C3 * (API / (temperature + 459.67)));
            }else{
                rso = O.C1 * gamma_100 * std::pow(Pb, O.C2) * (O.C3 * (API / (temperature + 459.67)));    
            }
        }
        return rso;
    }
    //Oil Bubble point pressure
    inline double Pb(const VBConstants_U30& U, const VBConstants_O30& O, double pressure, double temperature, double API, double Mg, double Psep, double Tsep, double Pb){
        //To be calculated with max Rso (meaning Rsb)
        double inside_var;
        double Pbubble;
        if(API < 30){
            inside_var = (Rso(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb) / (U.C1 * Mg * std::exp(U.C3 * (API / (temperature + 459.67))))); //removed Gamma_g function implementation
            Pbubble = std::pow(inside_var, (1 / U.C2));
        } else {
            inside_var = (Rso(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb) / (O.C1 * Mg * std::exp(O.C3 * (API / (temperature + 459.67))))); //removed Gamma_g function implementation
            Pbubble = std::pow(inside_var, (1 / O.C2));
        }
        return Pbubble;
    }

    //Oil FVF - Saturated Oil
    inline double bo_S(const VBConstants_U30& U, const VBConstants_O30& O, 
               double pressure, double temperature, double API, double Mg, double Psep, double Tsep, double Pb){
        
        double rso;
        double gamma_g;
        double bo;
        if(API < 30){
         rso = Rso(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb);
         gamma_g = Mg; //removed Gamma_g function implementation
         bo = 1 + (U.A1 * rso) + U.A2 * (temperature - 60) * (API / gamma_g) + U.A3 * (temperature - 60) * (API / gamma_g);
        } else{
         rso = Rso(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb);
         gamma_g = Mg; //removed Gamma_g function implementation
         bo = 1 + (O.A1 * rso) + O.A2 * (temperature - 60) * (API / gamma_g) + O.A3 * (temperature - 60) * (API / gamma_g);
        }
        return bo;
    }

    //Compressibility Undersaturated Oil
    inline double rsb_Pb(const VBConstants_U30& U, const VBConstants_O30& O, 
                         double Pb, double temperature, double API, double Mg, double Psep, double Tsep){
        double rsB = Rso(U, O, Pb, temperature, API, Mg, Psep, Tsep, Pb);

        return rsB;
    }
    inline double bo_B(const VBConstants_U30& U, const VBConstants_O30& O, 
                       double Pb, double temperature, double API, double Mg, double Psep, double Tsep){
        double boB = bo_S(U, O, Pb, temperature, API, Mg, Psep, Tsep, Pb);

        return boB;
        //Oil FVF for bubble point is the Bo value for saturated oil @ Pb
    }
    
    inline double co_U(const VBConstants_U30& U, const VBConstants_O30& O, 
                       double pressure, double temperature, double API, double Mg, 
                       double Psep, double Tsep, double rsB, double Pb){
        
        //If not provided with an Rsb (solution Rs @bubble point)
        if(rsB < 0){
            rsB = rsb_Pb(U, O, Pb, temperature, API, Mg, Psep, Tsep);
        }

        double gamma_g = Mg; //removed Gamma_g function implementation
        double coU = (-1433 + 5 * (rsB) + 17.2 * temperature - 1180 * gamma_g + 12.61 * API) / (pressure * std::pow(10, 5));

        return coU;
    }

    inline double bo_U(const VBConstants_U30& U, const VBConstants_O30& O, 
                       double pressure, double temperature, double API, double Mg, 
                       double Psep, double Tsep, double Pb){
        
        double rsB = rsb_Pb(U, O, Pb, temperature, API, Mg, Psep, Tsep);
        double boB = bo_S(U, O, Pb, temperature, API, Mg, Psep, Tsep, Pb);
        double coU = co_U(U, O, Pb, temperature, API, Mg, Psep, Tsep, rsB, Pb);

        double boU = boB * std::exp(coU * (Pb - pressure));

        return boU;
    
    }

    //Oil density Saturated - Undersaturated
    inline double rho_O(const VBConstants_U30& U, const VBConstants_O30& O, 
                       double pressure, double temperature, double API, double Mg, 
                       double Psep, double Tsep, double Pb, double rhoOSC_){
        
        double rhoO;
        if(pressure >= Pb){
            rhoO = rhoOSC_ / bo_U(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb);
        } else if (pressure < Pb){
            rhoO = rhoOSC_ / bo_S(U, O, pressure, temperature, API, Mg, Psep, Tsep, Pb);
        }

        return rhoO;
                       }

//-------------------4. Oil Viscosity Correlations -- Khan et al.-------------------

    /*
    Temperature   ---->> Rankine (R)
    Pressure      ---->>         (Psia)
    Viscocity     ---->>         (cp)
    */

    //!!!!!!!! ------------ gamma_o ------------ !!!!!!!!
        /* specific oil gravity for oil derived from API gamma_o = 141.5 / (131.5 + API)*/

    inline double theta_r(double temperature){

        double thetaR = (temperature + 459.67) / 459.67;

        return thetaR;

    }

    inline double muO_khPb(const VBConstants_U30& U, const VBConstants_O30& O, 
                         double Pb, double temperature, double API, double Mg, 
                         double Psep, double Tsep, double gamma_o){
        
        double gammaG = Mg; //removed Gamma_g function implementation
        double rsB = rsb_Pb(U, O, Pb, temperature, API, Mg, Psep, Tsep);
        double thetaR = theta_r(temperature);
        double inside_var = 1 - gamma_o;
        double muPb = (0.09 * std::sqrt(gammaG)) / (std::sqrt(rsB) * std::pow(thetaR, 4.5) * std::pow(inside_var, 3));

        return muPb;
        //bubble point viscocity
    }
    
    inline double muO_khU(const VBConstants_U30& U, const VBConstants_O30& O, double pressure, 
                         double Pb, double temperature, double API, double Mg, 
                         double Psep, double Tsep, double gamma_o){
        
        double muoB = muO_khPb(U, O, Pb, temperature, API, Mg, Psep, Tsep, gamma_o);
        double muoU = muoB * std::exp(9.6 * std::pow(10, -5) * (pressure - Pb));

        return muoU;
        //Undersaturated Oil viscocity
    }

    inline double muO_khS(const VBConstants_U30& U, const VBConstants_O30& O, double pressure, 
                         double Pb, double temperature, double API, double Mg, 
                         double Psep, double Tsep, double gamma_o){
        
        double muoB = muO_khPb(U, O, Pb, temperature, API, Mg, Psep, Tsep, gamma_o);
        double inside_var = pressure / Pb;
        double muoS = muoB * std::pow(inside_var, -0.14) * std::exp(-2.5 * std::pow(10, -4) * (pressure - Pb));

        return muoS;
        //Saturated Oil viscosity
    }
}


//Member initializer list for private members - Constructor
FluidProperties::FluidProperties(double defaultBo, double defaultBg, double defaultMuO,double defaultMuG, 
                                 double defaultMuW, double rhoORef, double rhoGRef, double rhoW, 
                                 double Pbubble, double API, double gasSG, double temperature, 
                                 double Pc, double Tc, double MgO, double MgG, double Psep, double Tsep, double gamma_o, double rhoO_SC)
    :defaultBo_(defaultBo), defaultBg_(defaultBg), defaultMuO_(defaultMuO), defaultMuG_(defaultMuG), 
    defaultMuW_(defaultMuW), rhoORef_(rhoORef), rhoGRef_(rhoGRef), rhoW_(rhoW), Pbubble_(Pbubble), 
    API_(API), gasSG_(gasSG), temperature_(temperature), Pc_(Pc), Tc_(Tc), MgO_(MgO), MgG_(MgG), Psep_(Psep), Tsep_(Tsep), gammaO_(gamma_o), rhoOSC_(rhoO_SC){

        //Convert Units
        double Psep_Psia = UnitConversion::Pa_to_Psia(Psep);
        double Tsep_F = UnitConversion::Kelvin_to_Fahrenheit(Tsep);
        double Tsep_R = UnitConversion::Kelvin_to_Rankine(Tsep);
        double Pb_Psia = UnitConversion::Pa_to_Psia(Pbubble);
        double temp_F = UnitConversion::Kelvin_to_Fahrenheit(temperature);
        double temp_R = UnitConversion::Kelvin_to_Rankine(temperature);
        double Tc_R = UnitConversion::Kelvin_to_Rankine(Tc);
        double Pc_Psia = UnitConversion::Pa_to_Psia(Pc);
        double MgO_lbmol = UnitConversion::gmol_to_lbmol(MgO);

        //Constants assignment
        Correlations::VBConstants_U30 vbU;
        Correlations::VBConstants_O30 vbO;
        Correlations::BWRConstants bwr;

        //Constant property defaults - lamda functions
        boFunc_= [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);

            //Ternary operator to determine Saturated VS Undersaturated
            return pPsia < Pb_Psia ? Correlations::bo_S(vbU, vbO, pPsia, temp_F, API, MgO_, Psep_Psia, Tsep_F, Pb_Psia) : 
                                     Correlations::bo_U(vbU, vbO, pPsia, temp_F, API, MgO_, Psep_Psia, Tsep_F, Pb_Psia);
        };
        bgFunc_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);
            double Bg = Correlations::computeBg_BWR8(pPsia, temp_R, Pc_Psia, Tc_R, bwr);
            return Bg;
        };
        
        muOFunc_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);

            //Ternary operator to determine Saturated VS Undersaturated
            return pPsia < Pb_Psia ? Correlations::muO_khS(vbU, vbO, pPsia, Pb_Psia, temp_R, API, MgO_, Psep_Psia, Tsep_R, gamma_o) 
                : Correlations::muO_khU(vbU, vbO, pPsia, Pb_Psia, temp_R, API, MgO_, Psep_Psia, Tsep_R, gamma_o);
        };

        muOBasicFunc_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);
            return Correlations::muO_khPb(vbU, vbO, pPsia, temp_R, API, MgO_, Psep_Psia, Tsep_R, gamma_o);
        };

        muGFunc_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);
            return Correlations::computeMuG_LGE(MgG_, pPsia, temp_R, Pc_Psia, Tc_R, bwr);
        };

        //Relative permeabilities (Corey-LIKE-curves)

            //Simplified Corey Parameters
                double krw0 = 0.3;
                double kro0 = 1;
                double krg0 = 0.6;
                double nw = 2;
                double no = 2;
                double ng = 2;
        
        krwFunc_= [&](double Sw){double krw = krw0 * std::pow(Sw, nw); return krw;};
        kroFunc_= [&](double So){double kro = kro0 * std::pow(So, no); return kro;};
        krgFunc_= [&](double Sg){double krg = krg0 * std::pow(Sg, ng); return krg;};

        //Rs correlation with bubblepoint
        rsFunc_= [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);
            return Correlations::Rso(vbU, vbO, pPsia, temp_F, API, MgO_, Psep_Psia, Tsep_F, Pb_Psia);
        };

        //Rho Oil density
        rhoO_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);
            return Correlations::rho_O(vbU, vbO, pPsia, Tc_R, API, MgO_, Psep_Psia, Tsep_R, Pb_Psia, rhoO_SC);
        };

        //Rho gas density
        rhoG_ = [=](double pressure){
            double pPsia = UnitConversion::Pa_to_Psia(pressure);

            return Correlations::computeRho_G(Pc_Psia, Tc_R, pPsia, temp_R, MgG_, bwr);
        };
}

//Formation Volume Factors
double FluidProperties::getBo(double pressure) const{
    return boFunc_(pressure);
}

double FluidProperties::getBg(double pressure) const{
    return bgFunc_(pressure);
}

//Viscosities
double FluidProperties::getMuO(double pressure) const{
    return muOFunc_(pressure);
}

double FluidProperties::getMuG(double pressure) const{
    return muGFunc_(pressure);
}

double FluidProperties::getMuW() const{
    return defaultMuW_;
}

double FluidProperties::getMuOBasic(double pressure) const{
    return muOBasicFunc_(pressure);
}

//Rs correlation
double FluidProperties::getRs(double pressure) const{
    return rsFunc_ ? rsFunc_(pressure) : 0.0;
}

//Densities
double FluidProperties::getRhoO(double pressure) const{
    
    return rhoO_(pressure); // Oil with dissolved gas
}

double FluidProperties::getRhoG(double pressure) const{
    return rhoG_(pressure);
}

double FluidProperties::getRhoW() const{
    return rhoW_;
}

//Relative Permeabilities
double FluidProperties::getKro(double So) const{
    return kroFunc_(So);
}
double FluidProperties::getKrw(double Sw) const{
    return krwFunc_(Sw);
}
double FluidProperties::getKrg(double Sg) const{
    return krgFunc_(Sg);
}

//Default Densities
double FluidProperties::getRhoOdef() const{
    return rhoORef_;
}

double FluidProperties::getRhoWdef() const{
    return rhoW_;
}

double FluidProperties::getRhoGdef() const{
    return rhoGRef_;
}

//Constants
double FluidProperties::getBubblePointPressure() const{
    return Pbubble_;
}