//
// Created by M. J. Steil on 2017.02.08.
//

#ifndef EOS_EOS_POLY_HPP
#define EOS_EOS_POLY_HPP

#include "eos.hpp"

#include <ctime>
#include <iomanip>

class eos_poly : public eos {
public:
    // Polytropic parameters derived from (kappa,gamma)
    double g, gm1, gm1inv, mBgm1invgK;
    double K, S,SnB;

    /**
     * Relativistic Polytrope constructor
     * @param kappa [1] dimensionless pressure coefficient
     * @param gamma [1] relativistic polyropic index gamma game in (1,2] for causal EoS at all h
     * @param mB_in = LmB_km [km] : mean baryon rest mass
     * @param h0_in = 0 [1] : h0
     * @param h1_in = 1.E-9 [1] : threshold for low h series expansions
     */
    eos_poly(double kappa, double gamma, double mB_in = LmB_km /*[km]*/, double h0_in = 0,  double h1_in = 1.E-9);

    /**
     * Setter method to set polytropic parameters
     * @param par = { kappa, gamma, mBin, h0, h1 }
     */
    void set_params(vector<double> par);

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

    // Save to file
    void save(string filename, int n, double h_min, double h_max);
};

#endif //EOS_EOS_POLY_HPP
