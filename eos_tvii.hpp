//
// Created by M. J. Steil on 2017.07.03.
//

#ifndef EOS_EOS_TVII_HPP
#define EOS_EOS_TVII_HPP

#include "eos.hpp"

#include <ctime>
#include <iomanip>
#include <complex>
#include <cmath>

using namespace std::literals::complex_literals;

class eos_tvii: public eos {
public:
    eos_tvii(double rhoC /*[km**-2]*/, double Z /*1*/, double mu /*1*/, double mBin=LmB_km /*[km]*/);

    // Main set method
    void set_params(double rhoC /*[km**-2]*/, double Z /*1*/, double mu /*1*/,double mBin=LmB_km /*[km]*/);
    void set_params(vector<double> par);

    // set with different parameters
    void set_RZ(double Rin, double Zin);
    void set_Z(double Zin);
    void set_hC(double hCin);

    // Global parameters of the Tolman VII configuration
    double rhoC, Z, mu;
    double M, R, RR, hC, PC, nbarC;

        // Auxiliary parameters
        double phi1, cosphi1, phi0, cosphi0 , C1, C2; // in reals
        double dPds0, drhods0, dnbards0, dhds0;

        complex<double> C3;
        complex<double> complex1 = 1.;

    // Auxiliary functions
    double Zofs(double s);

    complex<double> Wofs(double s);

    double phiofs(double s);
    double dphids(double s);

    // Thermodynamic quantities as function of the reduced interior radial variable s=r/R
    double Pofs(double s);
    double dPds(double s);

    double rhoofs(double s);
    double drhods(double s);

    double nbarofs(double s);
    double dnbards(double s);

    double hofs(double s);
    double dhds(double s);

    double csofs(double s);

    // sofh look up
    complex<double> Xofh(double h);
    double sofh(double h);

    // Metric potentials
    double lambdaofs(double s);
    double nuofs(double s);

    // Critical compactnesses
    double Zpinf(double muin = -1, int print = 0);
    double ZPrho(double muin = -1, double frac = 1., int print = 0);
    double Zcs(double muin = -1, double cs =1., int print = 0);

    // Export LORENE datafile
    void write_to_file(string filename);

    // Derived global parameters
    double MB(int print = 0);
    void star_info();

    // Thermodynamic quantities
    double P(double h);         // Pressure; [P]=km**-2
    double dlPdlh(double h);    // Logarithmic derivative of P; [dlPdlh]=1
    double dPdh(double h);      // Derivative of P; [dPdh]=[P]=km**-2

    double rho(double h);       // Energy density; [rho]=km**-2
    double dlrhodlh(double h);  // Logarithmic derivative of rho; [dlrhodlh]=1
    double drhodh(double h);    // Derivative of rho; [drhodh]=[rho]=km**-2

    double nbar(double h) ;     // Baryon number density; [nbar]=km**-3
    double dlnbardlh(double h); // Logarithmic derivative of nbar; [dlnbardlh]=1
    double dnbardh(double h);   // Derivative of rho; [dnbardh]=[nbar]=km**-3

    double h(double Qi, int i); // Logarithmic enthalpy h of the thermodynamic quantity Qi; i=0: P, i=1; rho, i=2: nbar
    double rhoofP(double P);    // Energy density; [rho]=km**-2 as function of pressure [km**-2]
};


#endif //EOS_EOS_TVII_TABLE
