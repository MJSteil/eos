//
// Created by M. J. Steil on 2017.02.08.
//

#ifndef EOS_EOS_IRF_HPP
#define EOS_EOS_IRF_HPP

#include "eos.hpp"

class eos_irf : public eos {
public:
    // IRF specific parameters
    double rhoC;

    /**
     * Incompressible relativistic fluid (IRF) constant density EoS constructor
     * @param rhoC_in [km**-2] : constant density
     * @param mB_in = LmB_km [km] : mean baryon rest mass
     */
    eos_irf(double rhoC_in /*[km**-2]*/, double mB_in = LmB_km /*[km]*/);

    /**
     * Setter method to set IRF parameters
     * @param par = { rhoC, mB, h0, h1 }
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

};


#endif //EOS_EOS_IRF_HPP
