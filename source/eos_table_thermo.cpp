//
// Created by M. J. Steil on 2017.07.04.
//

#include <gsl/gsl_interp.h>
#include "../eos_table.hpp"

// Thermodynamic quantities: Pressure P(h), energy density rho(h) and baryon number density nbar(h)
double eos_table::P(double h) {
    if(h<=h0) {
        return P0;
    }else{
        double lP_temp;
        double dlPdlh_temp;
        interpol_lP(h,lP_temp, dlPdlh_temp);

        return LrhoNuc_km*exp10(lP_temp); // [km**-2]
    }
}

double eos_table::rho(double h) {
    if(h<=h0){
        return rho0;
    }else{
        double lP_temp;
        double dlPdlh_temp;
        interpol_lP(h,lP_temp,dlPdlh_temp);
        double P_temp = exp10(lP_temp);

        return LrhoNuc_km*(P_temp/h * dlPdlh_temp - P_temp); // dP/dh = rho+P <=> rho = dP/dh-P = P/h*dlP/dlh-P, [km**-2]
    }
}

double eos_table::nbar(double h) {

    if(h<=h0){
        return nbar0;
    }else{

        double lP_temp;
        double dlPdlh_temp;
        interpol_lP(h,lP_temp,dlPdlh_temp);

        return (exp10(lP_temp) / h * dlPdlh_temp * exp(-h))*nbarScale; // exp(h) = (rho+P)/(mB nbar) <=> nbar = (rho+P)/mB*exp(-h) = P/h*dlP/dlh/mB*exp(-h), [km**-3]
    }

}

// Derivatives of the thermodynamic quantities P, rho and nbar
double eos_table::dPdh(double h) {
    if(h<=h0) {
        return P0 + rho0;
    }else{
        return P(h) + rho(h);
    }
}

double eos_table::drhodh(double h) {
    double dlPdlnbar_temp;
    interpol_lin_dlPdlnbar(h,dlPdlnbar_temp);

    if(h<=h0) {
        if(P0>0){
            return 1/dlPdlnbar_temp*dPdh(h)*(P0 + rho0)/P0;
        }else{
            return 0.;
        }
    }else{
        return 1/dlPdlnbar_temp*dPdh(h)*(P(h) + rho(h))/P(h);
    }
}

double eos_table::dnbardh(double h) {
    double dlPdlnbar_temp;
    interpol_lin_dlPdlnbar(h,dlPdlnbar_temp);

    if(h<=h0) {
        if(P0>0){
            return 1./dlPdlnbar_temp*nbar0*(1. + rho0/P0);
        }else{
            return 0.;
        }
    }else{
        return 1./dlPdlnbar_temp*nbar(h)*(1. + rho(h)/P(h));
    }

}


// Logarithmic derivatives of the thermodynamic quantities P, rho and nbar with dlQ/dlh=dQ/dh h/Q
double eos_table::dlPdlh(double h) {
    if(h<=h0) {
        return 0;
    }else{

        double lh_temp = log10(h);
        double lP_temp;
        double dlPdlh_temp;
        interpol_lP(h,lP_temp, dlPdlh_temp);

        return dlPdlh_temp;
    }
}

double eos_table::dlrhodlh(double h) {

    if(h<=h0){

        return dlnbardlh(h0);
    }else{

        return dlnbardlh(h)*(1. + P(h)/rho(h)); // dlrho/dlh = dlnbar/dlh*dlrho/dlnbar = dlnbar/dlh*(drho/dn n/rho) = dlnbar/dlh*(1+P/rho)
    }

}

double eos_table::dlnbardlh(double h) {
    if(h<=h0) {
        double dlPdlnbar_temp;
        interpol_lin_dlPdlnbar(h0,dlPdlnbar_temp);
        return 1./(dlPdlnbar_temp-1.);
    }else{

        double lP_temp;
        double dlPdlh_temp;
        double dlPdlnbar_temp;
        interpol_lP(h,lP_temp, dlPdlh_temp);
        interpol_lin_dlPdlnbar(h,dlPdlnbar_temp);

        return dlPdlh_temp/dlPdlnbar_temp; // dlnbar/dlh = dlnbar/dlP*dlP/dlh = dlP/dlh/(dlP/dlnbar)
    }
}


// Inverse function h(Qi,i)
double eos_table::h(double Qi, int i) {
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
                GSL_ERROR_VAL ("eos_table::h: called with invalid i. Choose i=0: P, i=1: rho, i=2: nbar.", GSL_EDOM,GSL_EDOM);
                hout = GSL_EDOM;
        }
    };

    gsl_function_pp<decltype(target_f)> star_target_fp(target_f);
    gsl_function *target_ff = static_cast<gsl_function*>(&star_target_fp);


    gsl_rootfinder star_target(target_ff, 60, 0,1E-16, {h0,0,hmax}, 0);
    hout = star_target.x_min;

    return hout;
}

// rho(P) using double logarithmic interpolation
double eos_table::rhoofP(double P) {
    double P_spline = log10(P/LrhoNuc_km);

    double rho_spline=gsl_spline_eval(rhoofPspline,P_spline,rhoofPsplineacc);

    return exp10(rho_spline)*cgcm3km2;
}