//
// Created by M. J. Steil on 2017.02.08.
//

#include "../eos_poly.hpp"

eos_poly::eos_poly(double kappa, double gamma, double mB_in, double h0_in,  double h1_in) {
    type    = "Relativistic Polytrope";
    set_params({kappa, gamma, mB_in, h0_in, h1_in});
}

void eos_poly::set_params(vector<double> par) {
    // Primary parameters
    g = par[1]; // Relativistic polyropic index
    K = (par[0])*LrhoNuc_km/(pow(LnNuc_km,g)); // Pressure coefficient with [K]=km**(3g-2)

    mB = par[2]; // Mean baryon rest mass; standard is mBL
    h0 = par[3]; // log(h)@P=0; standard is log(1)=0
    h1 = par[4]; // Threshold for series expansions log(h1); standard is 1E-9;

    // Auxiliary parameters
    gm1 = g-1.;
    gm1inv = 1./gm1;
    mBgm1invgK = mB*gm1/(K*g);

    // Scales
    SnB = pow(K,gm1inv);
    S = pow(K,gm1inv/2.);

    // Name
    stringstream name_tmp;
    name_tmp << "Polytropic EoS: (Kappa="<<par[0]<<"|Gamma="<<g<<"|mB="<<mB/cMeVkm<<" MeV)" ;
    name    = name_tmp.str();
}


// Thermodynamic quantities: Pressure P(h), energy density rho(h) and baryon number density nbar(h)
double eos_poly::P(double h) {
    if(h>h0){
        return K*pow(nbar(h),g);
    }else{
        return 0;
    }
}

double eos_poly::rho(double h) {
    if(h>h0){
        return gm1inv*P(h) + mB*nbar(h);
    }else{
        return 0;
    }
}

double eos_poly::nbar(double h) {
    if(h>h0){
        return pow(mBgm1invgK*(exp(h)-1),gm1inv);
    }else{
        return 0;
    }
}



// Derivatives of the thermodynamic quantities P, rho and nbar
double eos_poly::dPdh(double h) {
    if(h>h0){
        return dlPdlh(h)*P(h)/h;
    }else{
        return 0;
    }

}

double eos_poly::drhodh(double h) {
    if(h>h0){
        return dlrhodlh(h)*rho(h)/h;
    }else{
        return mB*mB/(2.*K)*(g==2.);
    }
}

double eos_poly::dnbardh(double h) {
    if(h>h0){
        return dlnbardlh(h)*nbar(h)/h;
    }else{
        return mB/(2.*K)*(g==2.);
    }
}

// Logarithmic derivatives of the thermodynamic quantities P, rho and nbar with dlQ/dlh=dQ/dh h/Q
double eos_poly::dlPdlh(double h) {
    // dlP/dh = dlP/dln*dln/dh = dlP/dn*n*dln/dh = g*dln/dh

    if(h>h0){
        if(h<h0+h1){
            return g*gm1inv*(1+h/2.+h*h/12.); // Series expansion around small h<h0+h1
        }else{
            return h*g*gm1inv/(1-exp(-h));
        }

    }else{
        return g*gm1inv;
    }

}

double eos_poly::dlrhodlh(double h) {
    // dlrho/dh = dlrho/dln*dln/dh = drho/dn*n*dln/dh = (1+P(h)/rho(h))*dln/dh

    if(h>h0){
        if(h<h0+h1){
            return (1+P(h)/rho(h))*gm1inv*(1+h/2.+h*h/12.); // Series expansion around small h<h0+h1
        }else{
            return (1+P(h)/rho(h))*h*gm1inv/(1-exp(-h));
        }

    }else{
        return gm1inv;
    }
}

double eos_poly::dlnbardlh(double h) {
    if(h>h0){
        if(h<h0+h1){
            return gm1inv*(1+h/2.+h*h/12.); // Series expansion around small h<h0+h1
        }else{
            return h*gm1inv/(1-exp(-h));
        }

    }else{
        return gm1inv;
    }
}

// Inverse function h(Qi,i)
double eos_poly::h(double Qi, int i) {
    double hout;
    switch (i)
    {
        case 0: // Ai=P
            hout =  log(1+pow(pow(Qi/K,1/g),gm1)/mBgm1invgK);
            break;
        case 1: // Ai=rho
            GSL_ERROR_VAL ("eos_poly::h: called with case i=1, which is not implemented. Choose i=0: P, i=2: nbar.", GSL_FAILURE,GSL_FAILURE);
            hout = GSL_FAILURE;
            break;
        case 2: // Ai=nbar
            hout =  log(1+pow(Qi,gm1)/mBgm1invgK);
            break;
        default:
            GSL_ERROR_VAL ("eos_poly::h: called with invalid i. Choose i=0: P, i=2: nbar.", GSL_EDOM,GSL_EDOM);
            hout = GSL_EDOM;
    }
    return hout;
}

double eos_poly::rhoofP(double P) {
    return mB*pow(P/K,1/g)+gm1inv*P;
}

// Save to file
void eos_poly::save(string filename, int n, double h_min, double h_max) {
    // File stream
    FILE *out;
    out=fopen(filename.c_str(),"w+");

    // Timestamp
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    stringstream time;
    time << put_time(&tm, "%Y.%m.%d - %H:%M");

    // lh-Grid
    double lhmin = log10(h_min);
    double lh_delta = (log10(h_max)-lhmin)/(n-1);

    function<double (int)> nbari = [&](int i)->double{
        double h = exp10(lhmin+lh_delta*i);
        return nbar(h)/cfm3km3;
    };

    function<double (int)> rhoi = [&](int i)->double{
        double h = exp10(lhmin+lh_delta*i);
        return rho(h)/cgcm3km2;
    };

    function<double (int)> Pi = [&](int i)->double{
        double h = exp10(lhmin+lh_delta*i);
        return P(h)/cdyncm2km2;
    };

    // Header
    fprintf(out,"# %s \n",time.str().c_str());
    fprintf(out,"# %s \n",name.c_str());
    fprintf(out,"# M. J. Steil \n");
    fprintf(out,"# Range: [%4E,%4E]fm**-3 \n",nbari(0),nbari(n-1));
    fprintf(out,"# \n");
    fprintf(out,"%d <-- Number of lines\n",n);
    fprintf(out,"#\n");
    fprintf(out,"#        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]\n");
    fprintf(out,"#\n");

    //Table
    for(int j=0;j<n;j++){
        if(j<9){
            fprintf(out,"    ");
        }else if(j<99){
            fprintf(out,"   ");
        }else {
            fprintf(out,"  ");
        }
        fprintf(out,"%d    %.15E    %.15E    %.15E\n",j+1,nbari(j),rhoi(j),Pi(j));
    }
    fclose(out);
}



