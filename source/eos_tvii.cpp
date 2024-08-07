//
// Created by M. J. Steil on 2017.07.03.
//

#include "../eos_tvii.hpp"

eos_tvii::eos_tvii(double rhoCin, double Zin, double muin, double mBin) {
    set_params(rhoCin,Zin,muin,mBin);
}

void eos_tvii::set_params(double rhoCin, double Zin, double muin, double mBin) {
    type = "Tolman VII solution/TVII EoS";

    mB   = mBin; // Mean baryon rest mass; standard is mBL
    h0  = 0; // log(h)@P=0; standard is log(1)=0
    h1  = 0; // Threshold for series expansions, NOT USED

    rhoC  = rhoCin;
    Z     = Zin;
    mu = muin;

    R = (sqrt(15/M_PI)*sqrt(Z))/(2.*sqrt(5*rhoC - 3*mu*rhoC));
    RR = R*R;
    M = Z*R;

    hC = log(1/cos(atan(sqrt((5*Z - 3*mu*Z)/(6*mu - 12*mu*Z))) + log((sqrt(2)*(-5 + 6*mu)*sqrt(Z) + 2*sqrt(3)*sqrt(mu*(-5 + 3*mu)*(-1 + 2*Z)))/(2*sqrt(3)*sqrt((5 - 3*mu)*mu) - 5*sqrt(2)*sqrt(Z)))/2.)/sqrt(1 + ((5 - 3*mu)*Z)/(6.*(mu - 2*mu*Z))));
    hmax=hC;
    
    phi1 = atan(sqrt(((5 - 3*mu)*Z)/(mu - 2*mu*Z))/sqrt(6.));
    cosphi1 = cos(phi1);
    phi0 = (2*atan(sqrt((5*Z - 3*mu*Z)/(6*mu - 12*mu*Z))) + log((sqrt(2)*(-5 + 6*mu)*sqrt(Z) + 2*sqrt(3)*sqrt(mu*(-5 + 3*mu)*(-1 + 2*Z)))/(2*sqrt(3)*sqrt((5 - 3*mu)*mu) - 5*sqrt(2)*sqrt(Z))))/2.;
    cosphi0 = cos(phi0);
    
    PC = (5*Z - sqrt(6)*sqrt(-(mu*(-5 + 3*mu)*Z))*tan(phi0))/(4.*(-5 + 3*mu)*M_PI*RR);
    nbarC = cosphi0/cosphi1*(PC+rhoC)/mB;
    
    C1 = (6*mu + 5*Z - 15*mu*Z)/(6.*mu);
    C2 = atan(sqrt(5 - 3*mu)/sqrt(mu*(-12 + 6/Z)));
    C3 = log((-5*sqrt(2) + 6*sqrt(2)*mu + 2*sqrt(3)*sqrt(complex1*mu*(-5 + 3*mu)*(2 - 1/Z)))/(6.*sqrt((5 - 3*mu)*mu)));

    dPds0 = (sqrt(mu*Z)*tan(phi0)*(-5*sqrt(6)*sqrt(5 - 3*mu)*Z + 3*(-5 + 3*mu)*sqrt(mu*Z)*tan(phi0)))/(2.*gsl_pow_2(5 - 3*mu)*M_PI*RR);
    drhods0 = (15*mu*Z)/(2.*(-5 + 3*mu)*M_PI*RR);
    dnbards0 = 1/(mB*exp(hC))*drhods0; // dn/ds= dn/drho drho/ds = 1/mu drho/ds = 1/(mB exp(h)) drho/ds

    dhds0 = (sqrt(6)*mu*sqrt((5 - 3*mu)/(mu*Z))*Z*tan(phi0))/(-5 + 3*mu);

    stringstream name_tmp;
    name_tmp << "TVII EoS: (rho="<<rhoC/cMeVfm3km2<<" MeVfm**-3 | hc= "<< hC << " | mB="<<mB/cMeVkm<<" MeV || Z="<< Z<<" | M="<< M/MSkm<<" MS | R="<< R << " km)" ;
    name    = name_tmp.str();
}

void eos_tvii::set_params(vector<double> par) {
    set_hC(par[0]);
}

// Auxiliary functions
double eos_tvii::Zofs(double s) {
    return (5.-10.*gsl_pow_2(s)*Z-3.*mu+6.*gsl_pow_4(s)*Z*mu)/(5.-3.*mu);
}

complex<double> eos_tvii::Wofs(double s) {
    return log((-5*sqrt(2) + 6*sqrt(2)*mu*gsl_pow_2(s) + 2*sqrt(3)*sqrt(complex1*(mu*(5 - 3*mu - 10*gsl_pow_2(s)*Z + 6*mu*gsl_pow_4(s)*Z))/Z))/(6.*sqrt((5 - 3*mu)*mu)));
}

double eos_tvii::phiofs(double s) {
    complex<double> phis = C2 + (C3-Wofs(s))*0.5;

    if(imag(phis)>0){
        char error_massage[100];
        sprintf(error_massage,"eos_tvii::phiof(s): imag(phiof(%.4E))=%.17E>0 !\n",s,imag(phis));
        GSL_ERROR_VAL (error_massage, GSL_FAILURE,1.);
    }
    return real(phis);
}

double eos_tvii::dphids(double s) {
    return - sqrt(6.)*s/sqrt((5.-3.*mu)*Zofs(s)/Z/mu); // also corrected in LORENE (PC 2018.06.26)
}

// Thermodynamic quantities as function of the reduced interior radial variable s=r/R
double eos_tvii::Pofs(double s) {
    return abs((s<=1)*((5. - 3.*mu*gsl_pow_2(s))*Z - sqrt(6)*tan(phiofs(s))*sqrt(-(mu*(-5. + 3.*mu)*Z*Zofs(s))))/(4.*(-5. + 3.*mu)*M_PI*RR));
}
double eos_tvii::dPds(double s) {
    return (s*sqrt(mu*Z)*tan(phiofs(s))*(3*(-5 + 3*mu)*sqrt(mu*Z)*tan(phiofs(s))*Zofs(s) + sqrt(6.)*(-5 + 6*mu*gsl_pow_2(s))*Z*sqrt((5 - 3*mu)*Zofs(s))))/(2.*gsl_pow_2(5 - 3*mu)*M_PI*RR*Zofs(s));
}


double eos_tvii::rhoofs(double s) {
    return abs(15./M_4PI/RR*Z*(mu*s*s-1)/(3*mu-5));
}
double eos_tvii::drhods(double s) {
    return 15.*s*Z*mu/(2*M_PI*RR*(-5.+3.*mu));
}


double eos_tvii::nbarofs(double s) {
    return abs(cos(phiofs(s))/cosphi1/mB*(Pofs(s)+rhoofs(s)));
}
double eos_tvii::dnbards(double s) {
    return ( cos(phiofs(s))*( drhods(s)+ dPds(s) ) - sin(phiofs(s))*(Pofs(s)+rhoofs(s))*dphids(s) )/mB/cosphi1;
}

double eos_tvii::hofs(double s) {
    return log(cosphi1/cos(phiofs(s)));
}
double eos_tvii::dhds(double s) {
    return tan(phiofs(s))*dphids(s);
}

double eos_tvii::csofs(double s) {
    if (s <= 1 && s >= 0) {
        return sqrt((tan(phiofs(s)) * (3 * (-5 + 3 * mu) * sqrt(mu * Z) * tan(phiofs(s)) * Zofs(s) +
                                       sqrt(6) * (-5 + 6 * mu * gsl_pow_2(s)) * Z * sqrt(-((-5 + 3 * mu) * Zofs(s))))) /
                    (15. * (-5 + 3 * mu) * sqrt(mu * Z) * Zofs(s)));
    } else {
        return 0;
    }
}

// sofh look up
complex<double> eos_tvii::Xofh(double h) {
    return exp(2.*acos(exp(-h)*cosphi1) - 2*C2 -C3);
}


double eos_tvii::sofh(double h){
    if(h==0){
        return 1;
    }else if(h>=hmax){
        return 0;
    }else{
        complex<double> sh = sqrt(60*mu + (18*sqrt(2)*sqrt(5 - 3*mu)*pow(mu,1.5))/Xofh(h) + (25*sqrt(2)*sqrt(mu)*Xofh(h))/sqrt(5 - 3*mu) - (6*sqrt(2)*sqrt(5 - 3*mu)*pow(mu,1.5)*Xofh(h))/Z)/(6.*sqrt(2)*mu);
        double sh_real = real(sh);
        double sh_imag = imag(sh);

        if(imag(sh)>0&&real(sh)>0){
            char error_massage[100];
            sprintf(error_massage,"eos_tvii::sof(h): imag(sof(%.4E))=%.4E>0 [%2E*real()] !\n",h,sh_imag,sh_imag/sh_real);
            //printf(ANSI_COLOR_YELLOW "%s" ANSI_COLOR_RESET,error_massage);
            GSL_ERROR_VAL (error_massage, GSL_FAILURE,GSL_FAILURE);
        }
        return sh_real;
    }
}

// Metric potentials
double eos_tvii::lambdaofs(double s) {
    return log(1./Zofs(s));
}

double eos_tvii::nuofs(double s) {
    return log(-((-6*mu + 5*(-1 + 3*mu)*Z)*gsl_pow_2(cos(phiofs(s))))/(6.*mu));
}

// Critical compactnesses
double eos_tvii::Zpinf(double muin, int print) {
    if(muin != -1){
        // Set a new mu if specified
        set_params(rhoC,Z,muin,mB);
    }

    // Finds phi(0)==M_PI/2.
    auto star_target_f = [&](double Z_in)->double{
        set_params(rhoC,Z_in,mu,mB);
        return  phiofs(0)-M_PI*0.5;
    };

    gsl_function_pp<decltype(star_target_f)> star_target_fp(star_target_f);
    gsl_function *star_target_ff = static_cast<gsl_function*>(&star_target_fp);


    gsl_rootfinder star_target(star_target_ff, 50, 1E-16,1E-16, {1E-8,0.,4/9.}, print);

    // Set found value in the star
    set_params(rhoC,star_target.x_min,mu,mB);

    return star_target.x_min;
}

double eos_tvii::ZPrho(double muin, double frac ,int print) {
    if(muin != -1){
        // Set a new mu if specified
        set_params(rhoC,Z,muin,mB);
    }

    // Zpinf of the config for upper bound for root finder
    double Zmax=Zpinf();

    // Finds P(0)/rho(0)==frac
    auto star_target_f = [&](double Z_in)->double{
        set_params(rhoC,Z_in,mu,mB);
        return  Pofs(0)/rhoofs(0)-frac;
    };

    gsl_function_pp<decltype(star_target_f)> star_target_fp(star_target_f);
    gsl_function *star_target_ff = static_cast<gsl_function*>(&star_target_fp);

    gsl_rootfinder star_target(star_target_ff, 50, 1E-16,1E-16, {1E-8,0.,Zmax}, print);

    // Set found value in the star
    set_params(rhoC,star_target.x_min,mu,mB);

    return star_target.x_min;
}

double eos_tvii::Zcs(double muin, double cs, int print) {
    if(muin != -1){
        // Set a new mu if specified
        set_params(rhoC,Z,muin,mB);
    }

    // Zpinf of the config for upper bound for root finder
    double Zmax=Zpinf();

    // Finds phi(0)==M_PI/2.
    auto star_target_f = [&](double Z_in)->double{
        set_params(rhoC,Z_in,mu,mB);
        return  csofs(0)-cs;
    };

    gsl_function_pp<decltype(star_target_f)> star_target_fp(star_target_f);
    gsl_function *star_target_ff = static_cast<gsl_function*>(&star_target_fp);

    gsl_rootfinder star_target(star_target_ff, 50, 1E-16,1E-16, {1E-8,0.,Zmax}, print);

    // Set found value in the star
    set_params(rhoC,star_target.x_min,mu,mB);

    return star_target.x_min;
}

// Derived global parameters
double eos_tvii::MB(int print) {

    auto comp_B_dBf = [&](double s)->double{
        return  s*s*nbarofs(s)*sqrt(exp(lambdaofs(s)));
    };
    gsl_function_pp<decltype(comp_B_dBf)> comp_B_dBfp(comp_B_dBf);
    gsl_function *comp_B_dB = static_cast<gsl_function*>(&comp_B_dBfp);

    // Integration Method
    gsl_integration comp_B_int(comp_B_dB,{0,1},1E-16,1E-16,"odeiv2",(size_t)5E4,6,print);

    if(print) printf("MB of tvii( mu=%.2f | rhoC=%.4E MeVfm**-3, Z=%.4E ) = %.16E MS",mu,rhoC/cMeVfm3km2,Z,M_4PI*R*RR*mB*comp_B_int.F/MSkm);

    return M_4PI*R*RR*mB*comp_B_int.F;
}

void eos_tvii::star_info() {
    printf("Tolman VII star:\n");
    printf("|-> P0=%.12E MeV/fm**3 <=> rho0=%.12E MeV/fm**3 <=> nbar0=%.8E fm**-3 <=> logh0=%.17g <=> cs=%.4E \n",PC/cMeVfm3km2,rhoC/cMeVfm3km2,nbarC/cfm3km3,hC,csofs(0));
    printf("|=> RS=%.17g km | MsG/Rs=%.17g | MsG=%.17g Msol | MsB=%.17g Msol\n",R,M/R,M/MSkm,MB(0)/MSkm);
    printf("\n");
}

void eos_tvii::set_RZ(double Rin, double Zin) {
    double rhoin = 15./M_4PI*Zin/(5 - 3*mu)/gsl_pow_2(Rin);

    set_params(rhoin,Zin,mu);
}

void eos_tvii::set_Z(double Zin) {
    set_params(rhoC,Zin,mu);
}

void eos_tvii::set_hC(double hCin) {
    // Zpinf of the config for upper bound for root finder
    double Zmax=Zpinf();
    double hmax = hC;

    // Finds phi(0)==M_PI/2.
    auto star_target_f = [&](double Z_in)->double{
        set_params(rhoC,Z_in,mu,mB);
        return  hCin-hC;
    };

    gsl_function_pp<decltype(star_target_f)> star_target_fp(star_target_f);
    gsl_function *star_target_ff = static_cast<gsl_function*>(&star_target_fp);

    gsl_rootfinder star_target(star_target_ff, 50, 1E-16,1E-16, {0.,0.,Zmax}, 0);

    // Set found value in the star
    set_params(rhoC,star_target.x_min,mu,mB);
    //printf("1-hCin/hc = %.8E\n",1-hCin/hC);

}

void eos_tvii::write_to_file(string filename) {
    FILE *out;
    out=fopen(filename.c_str(),"w+");

    fprintf(out,"7\tType of the EOS (cf. documentation of Eos::eos_from_file) \n");
    fprintf(out,"TVII%.0f\n",rhoC/cMeVfm3km2);
    fprintf(out,"%.16E\t// rhoC [km**-2]\n",rhoC);
    fprintf(out,"%.16E\t// Z=M/R [1]\n",Z);
    fprintf(out,"%.16E\t// mu [1]\n",mu);
    fprintf(out,"%.16E\t// mB [km]\n",mB);
    fprintf(out,"%.16E\t// ent0 [1]",h0);

    fclose(out);
}


