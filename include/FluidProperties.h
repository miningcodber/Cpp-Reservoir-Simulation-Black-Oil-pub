#pragma once
#include <vector>
#include <array>
#include <functional>

class FluidProperties{
    public:
        //Constructor
        FluidProperties(double defaultBo, double defaultBg,
                    double defaultMuO, double defaultMuG, double defaultMuW,
                    double rhoORef, double rhoGRef, double rhoW, double Pbubble, double API, double gasSG, double temperature,
                    double Pc, double Tc, double MgO, double MgG, double Psep, double Tsep, double gamma_o, double rhoO_SC);
        //Formation Volume Factors FVF
        double getBo(double pressure) const;
        double getBg(double pressure) const;

        //Viscosities
        double getMuO(double pressure) const;
        double getMuG(double pressure) const;
        double getMuW() const;
        double getMuOBasic(double pressure) const;

        //Relative permeabilities
        double getKro(double So) const;
        double getKrg(double Sg) const;
        double getKrw(double Sw) const;

        //Rs and densities
        double getRs(double pressure) const;
        double getRhoO(double pressure) const;
        double getRhoG(double pressure) const;
        double getRhoW() const;

        //Default Densities
        double getRhoOdef() const;
        double getRhoWdef() const;
        double getRhoGdef() const;

        //Constants
        double getBubblePointPressure() const;

        
    private:
        //Default constants
        double defaultBo_, defaultBg_;
        double defaultMuO_, defaultMuG_, defaultMuW_;
        double rhoORef_, rhoGRef_, rhoW_;
        double Pbubble_;
        double API_;
        double gasSG_;
        double temperature_;
        double Pc_;
        double Tc_;
        double Mw_;
        double Psep_;
        double Tsep_;
        double gammaO_;
        double rhoOSC_;
        double MgO_;
        double MgG_;


        //Optional functional override functions
        std::function<double(double)> boFunc_, bgFunc_;
        std::function<double(double)> muGFunc_, muOFunc_, muOBasicFunc_;
        std::function<double(double)> kroFunc_, krgFunc_, krwFunc_;
        std::function<double(double)> rsFunc_;
        std::function<double(double)> rhoO_;
        std::function<double(double)> rhoG_;

};