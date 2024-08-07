//
// Created by M. J. Steil on 2017.02.09.
//

#ifndef EOS_EOS_TABLE_HPP
#define EOS_EOS_TABLE_HPP

#include "eos.hpp"

#include <gsl/gsl_spline.h>

class eos_table : public eos {
public:
    eos_table();
    eos_table(string filename_in, string type_in, vector<double> params = {}, int print = 0);

    // Read files
    void read_cgs(string filename_in, int print = 0);
        // EoS Raw data
        string filename, info;
        size_t table_length;
        vector<int> no_data;
        vector<double> nbar_fm3_data, rho_cgs_data, P_cgs_data;

    // Generate EoS data
    void generate_data( int print = 0);
        // EoS data
        double P0, rho0, nbar0;
        double mu0scale;
        double nbarScale;

        vector<double> h_data;
        vector<double> lP_data, dlP_data;
        vector<double> dlPdlnbar_data;

        // EoS Auxiliary data
        double lhth, hth ,mu0;
        double lh1=-1.E200;
        double lh2=-1.E200;
        int h1h2muManual = 0;
        vector<double> lh_data;
        vector<double> lnbar_data, lrho_data; // references

    // Generate interpolations
    void generate_int();
        // Interpolation
        gsl_interp_accel *lpacc;
        void interpol_herm_lP(double lh,double &lP, double &dlP);

        gsl_interp_accel *lp_low_acc;
        void interpol_lP_low(double h,double &lP, double &dlP);

        void interpol_lP(double h,double &lP, double &dlP);

        gsl_interp_accel *dlPdlnbaracc;
        gsl_spline *dlPdlnbarspline;
        void interpol_lin_dlPdlnbar(double h,double &dlPdlnbar);

        // rho(P) direct Log-Log-Spline
        gsl_interp_accel *rhoofPsplineacc;
        gsl_spline *rhoofPspline;

    // Inherited setter: Has no function for eos_tab
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

    void plot_interpolation(int plot);
    void read_cgs_old(string filename_in, double mu0 =0., int print=0);
};

#endif //EOS_EOS_TABLE_HPP
