//
// Created by M. J. Steil on 2017.02.08.
//

#include "../eos_irf.hpp"

eos_irf::eos_irf(double rhoC_in /*[km**-2]*/, double mB_in /*= LmB_km [km]*/) {
    type = "Incompressible relativistic fluid (IRF) constant density EoS";

    set_params({rhoC_in,mB_in,0,0});
}

void eos_irf::set_params(vector<double> par) {
    rhoC = par[0]; // Constant energy density of the IRF EoS

    mB  = par[1]; // Mean baryon rest mass; standard is mBL
    h0  = par[2]; // log(h)@P=0; standard is log(1)=0
    h1  = par[3]; // Threshold for series expansions, NOT USED

    stringstream name_tmp;
    name_tmp << "ISS - Constant density EoS: (rho="<<rhoC/cMeVfm3km2<<" MeVfm**-3 | mB="<<mB/cMeVkm<<" MeV)" ;
    name    = name_tmp.str();
}


// Thermodynamic quantities: Pressure P(h), energy density rho(h) and baryon number density nbar(h)
double eos_irf::P(double h) {
    if(h>h0){
        return rhoC*(exp(h)-1);
    }else{
        return 0;
    }
}

double eos_irf::rho(double h) {
    return rhoC;
}


double eos_irf::nbar(double h) {
    return rhoC/mB;
}

// Derivatives of the thermodynamic quantities P, rho and nbar
double eos_irf::dPdh(double h) {
    if(h>h0){
        return rhoC*exp(h);
    }else{
        return 0;
    }
}

double eos_irf::drhodh(double h) {
    return 0;
}

double eos_irf::dnbardh(double h) {
    return 0;
}

// Logarithmic derivatives of the thermodynamic quantities P, rho and nbar with dlQ/dlh=dQ/dh h/Q
double eos_irf::dlPdlh(double h) {
    if(h>h0){
        return exp(h)*h/(exp(h)-1.);  // dlP/dlh = dP/dh*P/h
    }else{
        return 0;
    }
}

double eos_irf::dlrhodlh(double h) {
    return 0;
}

double eos_irf::dlnbardlh(double h) {
    return 0;
}


// Inverse function h(Qi,i)=h(P) for i=0
double eos_irf::h(double Qi, int i) {
   if(i==0){
       return log(Qi/rhoC+1);
   }else{
       GSL_ERROR_VAL ("eos_irf::h: called with invalid i. Choose i=0: P.", GSL_EDOM,GSL_EDOM);
   }
}

// rho(P)
double eos_irf::rhoofP(double P) {
    return rhoC;
}

