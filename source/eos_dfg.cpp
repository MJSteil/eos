//
// Created by M. J. Steil on 2017.08.16.
//

#include "../eos_dfg.hpp"

eos_dfg::eos_dfg() {
    type    = "Degenerate Fermi gas EoS";
    set_params({mn*cMeVkm});
}

eos_dfg::eos_dfg(double mB_in) {
    type    = "Degenerate Fermi gas EoS";
    set_params({mB_in});
}

void eos_dfg::set_params(vector<double> par) {
    // Primary parameters
    mB=par[0];
    mBhbar3=gsl_pow_3(mB/(hbar*G/c0_3*1E-6))/M_PI_SQARE; // [km**-3]
    mBhbar4=mB*mBhbar3; // [km**-2]

    h0=0;

    // Name
    stringstream name_tmp;
    name_tmp << "Degenerate Fermi gas EoS: (mB="<<mB/cMeVkm<<" MeV)" ;
    name    = name_tmp.str();
}


// Thermodynamic quantities: Pressure P(h), energy density rho(h) and baryon number density nbar(h)
double eos_dfg::P(double h) {
    if(h>h0){
        double exph = exp(h);
        double exp2h = exp(2.*h);
        return mBhbar4/24.*( exph*sqrt(exp2h-1.)*(2.*exp2h-5.) + 3*log(exph+sqrt(exp2h-1)) );
    }else{
        return 0;
    }
}

double eos_dfg::rho(double h) {
    if(h>h0){
        double exph = exp(h);
        double exp2h = exp(2.*h);
        return mBhbar4/8.*( exph*(2.*exp2h-1.)*sqrt(exp2h-1.) - log(exph+sqrt(exp2h-1.)) );
    }else{
        return 0;
    }
}

double eos_dfg::nbar(double h) {
    if(h>h0){
        return mBhbar3/3.*pow(exp(2.*h)-1,1.5);
    }else{
        return 0;
    }
}



// Derivatives of the thermodynamic quantities P, rho and nbar
double eos_dfg::dPdh(double h) {
    if(h>h0){
        return mBhbar4/3.*pow(exp(2.*h)-1,1.5)*exp(h);
    }else{
        return 0;
    }

}

double eos_dfg::drhodh(double h) {
    if(h>h0){
        return mBhbar4*sqrt(exp(2.*h)-1.)*exp(3.*h);
    }else{
        return 0;
    }
}

double eos_dfg::dnbardh(double h) {
    if(h>h0){
        return mBhbar3*exp(2.*h)*sqrt(exp(2.*h)-1.);
    }else{
        return 0;
    }
}

// Logarithmic derivatives of the thermodynamic quantities P, rho and nbar with dlQ/dlh=dQ/dh h/Q
double eos_dfg::dlPdlh(double h) {
    if(h>h0){
        return dPdh(h)/P(h)*h;
    }else{
        return 0;
    }

}

double eos_dfg::dlrhodlh(double h) {
    if(h>h0){
        return drhodh(h)/rho(h)*h;
    }else{
        return 0;
    }
}

double eos_dfg::dlnbardlh(double h) {
    if(h>h0){
        return dnbardh(h)/nbar(h)*h;
    }else{
        return 0;
    }
}

// Inverse function h(Qi,i)
double eos_dfg::h(double Qi, int i) {
    double hout;

    auto target_f = [&](double hC_in)->double{
        switch (i)
        {
            case 0: // Ai=P
                return  Qi-P(hC_in);
            case 1: // Ai=rho
                return  Qi-rho(hC_in);
            case 2: // Ai=nbar
                return  Qi-nbar(hC_in);
            default:
                GSL_ERROR_VAL ("eos_dfg::h: called with invalid i. Choose i=0: P, i=1: rho, i=2: nbar.", GSL_EDOM,GSL_EDOM);
                hout = GSL_EDOM;
        }
    };

    gsl_function_pp<decltype(target_f)> star_target_fp(target_f);
    gsl_function *target_ff = static_cast<gsl_function*>(&star_target_fp);


    gsl_rootfinder star_target(target_ff, 60, 0,1E-16, {h0,0,1E2}, 0);
    hout = star_target.x_min;

    return hout;
}

double eos_dfg::rhoofP(double P) {
    return rho(h(P,0));
}



