//
// Created by M. J. Steil on 2017.07.04.
//

#include "../eos_tvii.hpp"

// Thermodynamic quantities: Pressure P(h), energy density rho(h) and baryon number density nbar(h)
double eos_tvii::P(double h) {
    double sh = sofh(h);
    return Pofs(sh);
}

double eos_tvii::rho(double h) {
    double sh = sofh(h);
    return rhoofs(sh);
}

double eos_tvii::nbar(double h) {
    double sh = sofh(h);
    return nbarofs(sh);
}

// Derivatives of the thermodynamic quantities P, rho and nbar
double eos_tvii::dPdh(double h) {
    double sh = sofh(h);
    if(sh!=0){
        return dPds(sh)/dhds(sh);
    }else{
        return dPds0/dhds0;
    }

}

double eos_tvii::drhodh(double h) {
    double sh = sofh(h);
    if(sh!=0){
        return drhods(sh)/dhds(sh);
    }else{
        return drhods0/dhds0;
    }
}

double eos_tvii::dnbardh(double h) {
    double sh = sofh(h);
    if(sh!=0){
        return dnbards(sh)/dhds(sh);
    }else{
        return dnbards0/dhds0;
    }
}

// Logarithmic derivatives of the thermodynamic quantities P, rho and nbar with dlQ/dlh=dQ/dh h/Q
double eos_tvii::dlPdlh(double h) {
    return dPdh(h)*h/P(h);
}

double eos_tvii::dlrhodlh(double h) {
    return drhodh(h)*h/rho(h);
}

double eos_tvii::dlnbardlh(double h) {
    return dnbardh(h)*h/nbar(h);
}


// Inverse function h(Qi,i)
double eos_tvii::h(double Qi, int i) {
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
                GSL_ERROR_VAL ("eos_poly::h: called with invalid i. Choose i=0: P, i=2: nbar.", GSL_EDOM,GSL_EDOM);
                hout = GSL_EDOM;
        }
    };

    gsl_function_pp<decltype(target_f)> star_target_fp(target_f);
    gsl_function *target_ff = static_cast<gsl_function*>(&star_target_fp);


    gsl_rootfinder star_target(target_ff, 60, 0,1E-18, {h0,0,hmax}, 0);
    hout = star_target.x_min;

    return hout;
}

double eos_tvii::rhoofP(double P) {
    GSL_ERROR_VAL ("eos_tvii::rhoofP: Not implemented.", GSL_FAILURE, GSL_FAILURE);
}
