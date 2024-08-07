/***
 * eos
 * Created by M. J. Steil on 2017.01.19.
 */

#ifndef EOS_HPP
#define EOS_HPP

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>

#include "../units/units.hpp"
#include "../gsl_wrapper/gsl_wrapper.hpp"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;
using namespace units;
using namespace gsl_wrapper;

class eos {
public:
    string name;
    string type = "EoS parent class";

    double mB   = 0.0f; // Mean baryon rest mass; standard is mBL
    double h0  = 0.0f; // log(h)@P=0; standard is log(1)=0
    double hmax  = 0.0f; // Maximal value for h
    double h1  = 0.0f; // Threshold for series expansions

    eos();

    virtual void set_params(vector<double> par) = 0; // Generic setter method to adjust various EoS parameters on the fly

    virtual double P(double h)      = 0; // Pressure; [P]=km**-2
    virtual double dlPdlh(double h) = 0; // Logarithmic derivative of P; [dlPdlh]=1
    virtual double dPdh(double h)   = 0; // Derivative of P; [dPdh]=[P]=km**-2

    virtual double rho(double h)        = 0; // Energy density; [rho]=km**-2
    virtual double dlrhodlh(double h)   = 0; // Logarithmic derivative of rho; [dlrhodlh]=1
    virtual double drhodh(double h)     = 0; // Derivative of rho; [drhodh]=[rho]=km**-2

    virtual double nbar(double h)       = 0; // Baryon number density; [nbar]=km**-3
    virtual double dlnbardlh(double h)  = 0; // Logarithmic derivative of nbar; [dlnbardlh]=1
    virtual double dnbardh(double h)    = 0; // Derivative of rho; [dnbardh]=[nbar]=km**-3

    virtual double h(double Qi, int i)  = 0; // Logarithmic enthalpy h of the thermodynamic quantity Qi; i=0: P, i=1; rho, i=2: nbar
    virtual double rhoofP(double P)     = 0; // Energy density; [rho]=km**-2 as function of pressure [km**-2]

    // Non virtual eos base class methods
    double cs(double h, double hthM =1E1);
    void compare(eos &eos_in, double h); // Compare eos to eos_in
    void info();
    void nbar_sat(double h);
};

#endif //EOS_HPP
