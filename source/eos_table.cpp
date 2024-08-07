//
// Created by M. J. Steil on 2017.02.09.
//

#include "../eos_table.hpp"
// Helper
std::istream& safeGetline2(std::istream& is, std::string& t) {
    // Source: http://stackoverflow.com/a/6089413, 170209 15:08
    t.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

// Setup
eos_table::eos_table() {

}

eos_table::eos_table(string filename_in, string type_in, vector<double> params, int print) {

    if(params.empty()-1){
        mu0 = params[0];
        lh2 = params[1];
        lh1 = params[2];
        h1h2muManual = 1;
    }

    if(type_in=="cgs"){
        read_cgs(filename_in,print);
        generate_data(print);

        generate_int();

        if(print) eos::info();
    }else if(type_in=="cgs_old"){
        read_cgs_old(filename_in,0,print);
        generate_int();

        if(print) eos::info();
    }else{
        GSL_ERROR_VAL ("eos_table: Invalid type_in! At the moment 'cgs' only.", GSL_FAILURE,);
    }
}

void eos_table::set_params(vector<double> par) {

}

void eos_table::read_cgs(string filename_in, int print) {
    type    = "CGS tabulated EoS";
    filename = filename_in;

    //region Open file: ifstream fich
    ifstream fich(filename.c_str()) ;
    if (!fich) {
        cerr << "eos_table: " << "Can not open " << filename.c_str() <<"!"<< endl ;
        cerr << "Aborting..." << endl ;
        abort() ;
    }
    //endregion

    //region Process Header
    /* Header format example:
    #  2017.04.24 - 14:10
    #  Hempel-Schaffner-Bielich/DD2 [HSDD2]
    #  M. Hempel and J. Schaffner-Bielich, Nucl. Phys. A 837 (2010) 210
    #  http://compose.obspm.fr/spip.php?article21 - 2017.04.24 - 14:00 (new link http://compose.obspm.fr/eos/18/)
    #  { 1.000E-9,  1.660E6,  9.490E22} to { 1.000E1,  1.550E17,  1.280E38}
    247     0     -1.7    -4   <-- Number of lines (cout out data points 102-105), mu0, hint_O2_th, hint_O1_th
    #
    #        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
    #
    */
    fich.ignore(1000, '\n') ; // Jump over first line
    fich.ignore(1);
    safeGetline2(fich, info);// Get info (header line 2)

    for (int i=0; i<3; i++) {fich.ignore(1000, '\n');}; // Skip to line 6

    fich >> table_length ;  // Get table_length
    if(h1h2muManual!=1) {
        fich >> mu0;   // mu0
        fich >> lh2;    // lh1 threshold
        fich >> lh1;    // lh2 threshold
    }
    fich.ignore(1000, '\n') ;

    for (int i=0; i<3; i++) {fich.ignore(1000, '\n');}; // Skip column headers
    //endregion

    mB  = LmB_km;  // Mean baryon rest mass; standard is mBL, not used by this class directly

    // Setup raw data vectors
    no_data         = vector<int>(table_length); 
    nbar_fm3_data   = vector<double>(table_length);
    P_cgs_data      = vector<double>(table_length);
    rho_cgs_data    = vector<double>(table_length);

    // Read in data and close stream
    if(print) printf("-: no | nbar [fm**-3], rho [g cm**-3], P [dyne cm**-2] \n");
    for (int i=0; i<table_length; i++) {
        fich >> no_data[i] ;
        fich >> nbar_fm3_data[i];
        fich >> rho_cgs_data[i]  ;
        fich >> P_cgs_data[i];
        fich.ignore(1000, '\n') ;


        if(print) printf("-: %d | %.4E, %.4E, %.4E \n",no_data[i],nbar_fm3_data[i],rho_cgs_data[i] ,P_cgs_data[i]);
    }
    fich.close();
    
    // Generate table name string
    stringstream name_tmp;
    name_tmp << "Tabulated EoS (CGS):" << info.c_str() <<" - (from: " << filename.c_str() <<", "<< table_length <<" points)";
    name    = name_tmp.str();
}

void eos_table::generate_data(int print) {
    //region Process in raw data
    double nbar_tmp, rho_tmp, P_tmp, h_tmp;

    // Setup main data containers
    h_data  = vector<double>(table_length);
    lP_data     = vector<double>(table_length);
    dlP_data    = vector<double>(table_length);

    lnbar_data  = vector<double>(table_length);
    lrho_data   = vector<double>(table_length);
    
    // Write into main data containers
    if(print)printf("-: i | lP_data[i], dlP_data[i], lrho_data[i], lnbar_data[i] | h_data[i], h_data[i]-h_data[i-1] \n");
    for (int i=0; i<table_length; i++) {

        // Interpolation is done for P and rho in CGS and for nbar in fm**-3 to be consistent with LORENE eos_tabul
        P_tmp       = P_cgs_data[i]/c0c0_cgs;   // [P_tmp] = gcm**-3
        nbar_tmp    = nbar_fm3_data[i];         // [nbar_tmp] = fm**-3
        rho_tmp     = rho_cgs_data[i];          // [rho_tmp] = gcm**-3

        // mu0/h0-Setup
        if(i==0) {
            if (mu0 == 0) {
                mu0scale = (P_tmp + rho_tmp) / (nbar_tmp); // Natural choice for h(P(0))=0;
            } else {
                mu0scale = mu0;
            }
            nbarScale = 1 / mu0scale * LrhoNuc_cgs * cfm3km3;
        }
        
        h_tmp = log((P_tmp+rho_tmp)/(nbar_tmp*mu0scale));

        h_data[i] = h_tmp;
        lP_data[i] = log10(P_tmp/LrhoNuc_cgs);
        dlP_data[i] = (rho_tmp+P_tmp)/P_tmp*h_tmp;

        lnbar_data[i] = log10(nbar_tmp);
        lrho_data[i] = log10(rho_tmp);
        
        if(i>0){
            if(print)printf("-: %d | %.4E, %.4E, %.4E, %.4E | %.4E, %d \n",i,lP_data[i],dlP_data[i],lrho_data[i],lnbar_data[i],h_data[i],(h_data[i]>h_data[i-1]));
        }
    }

    h0 = h_data[0];
    P0 = exp10(lP_data[0])*LrhoNuc_km;
    rho0 = exp10(lrho_data[0])*cgcm3km2;
    nbar0 = exp10(lnbar_data[0])*cfm3km3;
    hmax = h_data[table_length-1];
    //endregion

    //region Compute dlPdlnbar_data using cubic polynomials and setup dlPdlnbarspline
    dlPdlnbar_data = vector<double> (table_length);
    vector<double> p3(3);
    vector<double> n3(3);

    // special case: i=0 - Surface
    p3 = { lP_data[0], lP_data[1], lP_data[2] };
    n3 = { lnbar_data[0], lnbar_data[1], lnbar_data[2] };


    dlPdlnbar_data[0] = p3[0]*(2*n3[0]-n3[1]-n3[2])/(n3[0]-n3[1])/(n3[0]-n3[2]) +
                        p3[1]*(n3[0]-n3[2])/(n3[1]-n3[0])/(n3[1]-n3[2]) +
                        p3[2]*(n3[0]-n3[1])/(n3[2]-n3[0])/(n3[2]-n3[1]) ;

    // Main loop 1<=i<table_length-1
    for(int i=1;i<table_length-1;i++) {
        p3 = { lP_data[i-1], lP_data[i], lP_data[i+1] };
        n3 = { lnbar_data[i-1], lnbar_data[i], lnbar_data[i+1] };


        dlPdlnbar_data[i] = p3[0]*(n3[1]-n3[2])/(n3[0]-n3[1])/(n3[0]-n3[2]) +
                            p3[1]*(2*n3[1]-n3[0]-n3[2])/(n3[1]-n3[0])/(n3[1]-n3[2]) +
                            p3[2]*(n3[1]-n3[0])/(n3[2]-n3[0])/(n3[2]-n3[1]) ;
    }

    // special case: i=table_length-1 - Nucleus
    p3 = { lP_data[table_length-3], lP_data[table_length-2], lP_data[table_length-1] };
    n3 = { lnbar_data[table_length-3], lnbar_data[table_length-2], lnbar_data[table_length-1] };


    dlPdlnbar_data[table_length-1] = p3[0]*(n3[2]-n3[1])/(n3[0]-n3[1])/(n3[0]-n3[2]) +
                                     p3[1]*(n3[2]-n3[0])/(n3[1]-n3[0])/(n3[1]-n3[2]) +
                                     p3[2]*(2*n3[2]-n3[0]-n3[1])/(n3[2]-n3[0])/(n3[2]-n3[1]);

    //endregion
}

void eos_table::generate_int() {
    hth = (h_data[1]);
    lhth = log10(hth);

    vector<double> lh_in(table_length);
    for(int i=1; i<table_length; i++){
        lh_in[i]= log10(h_data[i]);
    }
    lh_in[0]=lh_in[1];
    lh_data.assign(lh_in.begin(),lh_in.begin()+table_length);

    // interpol_herm_lP accelerator
    lpacc = gsl_interp_accel_alloc();

    // interpol_lP_low accelerator
    lp_low_acc = gsl_interp_accel_alloc();

    // Linear spline setup
    dlPdlnbaracc    = gsl_interp_accel_alloc ();
    dlPdlnbarspline = gsl_spline_alloc (gsl_interp_linear, table_length);
    gsl_spline_init (dlPdlnbarspline, h_data.data(), dlPdlnbar_data.data(), table_length);

    // rho(P) Log-Log-Spline setup
    rhoofPsplineacc    = gsl_interp_accel_alloc ();
    rhoofPspline = gsl_spline_alloc (gsl_interp_cspline, table_length);
    gsl_spline_init (rhoofPspline, lP_data.data(), lrho_data.data(), table_length);

}

// Interpolations
void eos_table::interpol_herm_lP(double lh,double &lP, double &dlP) {

    // Interpolation index lookup - gsl_interp_accel_find
    size_t i = gsl_interp_accel_find (lpacc, lh_data.data(), table_length, lh);
    size_t ip1 = i+1;

    double dx = lh_data[ip1]-lh_data[i];
    double dxdx = dx*dx;

    double x1   = lh_data[i];
    double y1   = lP_data[i];
    double y2   = lP_data[ip1];
    double dy1  = dlP_data[i];
    double dy2  = dlP_data[ip1];

    double u = (lh - lh_data[i]) / dx ;
    double u2 = u*u ;
    double u3 = u2*u ;


    if(lh>lh2&&lh>lh1){

        // Cubic hermite interpolation polynomial of O[lh**3]
        lP =   y1 * ( 2.*u3 - 3.*u2 + 1.)
               + y2 * ( 3.*u2 - 2.*u3)
               + dy1 * dx * ( u3 - 2.*u2 + u )
               - dy2 * dx * ( u2 - u3 ) ;

        dlP =   6. * ( y1 - y2 ) * ( u2 - u ) / dx
                + dy1* ( 3.*u2 - 4.*u + 1.)
                + dy2 * ( 3.*u2 - 2.*u );

    }else if(lh>lh1){

        // Quadratic interpolation polynomial of O[lh**2]
        lP = y1 + u*(dx*dy1 + u*(-(dx*dy1) - y1 + y2));
        dlP = dy1 + (-dy1 + dy2)*u;

    }else{
        // Linear interpolation polynomial of O[lh**2]
        lP = y1 + u*(-y1 + y2);
        dlP =dy1 + (-dy1 + dy2)*u;

    }
}

void eos_table::interpol_lP_low(double h,double &lP, double &dlP) {

    // Interpolation index lookup - gsl_interp_accel_find
    size_t i = gsl_interp_accel_find (lp_low_acc, h_data.data(), table_length, h);
    size_t ip1 = i+1;

    double dx = h_data[ip1]-h_data[i];
    double dxdx = dx*dx;

    double x1   = h_data[i];
    double y1   = lP_data[i];
    double y2   = lP_data[ip1];
    double dy1  = dlP_data[i];
    double dy2  = dlP_data[ip1];

    double u = (h - h_data[i]) / dx ;
    double u2 = u*u ;
    double u3 = u2*u ;


    // Linear interpolation polynomial of O[lh**2]
    lP = y1 + u*(-y1 + y2);
    dlP =dy1 + (-dy1 + dy2)*u;

/*    // Quadratic interpolation polynomial of O[h**2]
    lP = y1 + u*(dx*dy1 + u*(-(dx*dy1) - y1 + y2));
    double dP_tmp = dy1 + (-dy1 + dy2)*u;

    if(h>0){
        dlP =dP_tmp;
    }else{
        dlP =0;
    }*/

/*    // Cubic hermite interpolation polynomial of O[lh**3]
    lP =   lP_data[i] * ( 2.*u3 - 3.*u2 + 1.)
           + lP_data[ip1] * ( 3.*u2 - 2.*u3)
           + dlP_data[i] * dx * ( u3 - 2.*u2 + u )
           - dlP_data[ip1] * dx * ( u2 - u3 ) ;

    dlP =   6. * ( lP_data[i] - lP_data[ip1] ) * ( u2 - u ) / dx
            + dlP_data[i] * ( 3.*u2 - 4.*u + 1.)
            + dlP_data[ip1] * ( 3.*u2 - 2.*u );


    // Linear interpolation polynomial of O[lh**2]
    lP = y1 + u*(-y1 + y2);
    dlP =dy1 + (-dy1 + dy2)*u;
    //cout << "1" <<endl;*/

}

void eos_table::interpol_lP(double h,double &lP, double &dlP) {
    //region bounds error
    if(h>hmax){
        GSL_ERROR_VAL ("interpol_lin_dlPdlnbar: lh>lhmax", GSL_EDOM,);
    }
    //endregion

    if(h<=hth){
        interpol_lP_low(h,lP, dlP);
    }else{
        double lh_temp = log10(h);
        interpol_herm_lP(lh_temp,lP, dlP);
    }

}

void eos_table::interpol_lin_dlPdlnbar(double h, double &dlPdlnbar) {

    //region lh>lhmax error
//    if(h>hmax){
//        GSL_ERROR_VAL ("interpol_lin_dlPdlnbar: lh>lhmax", GSL_EDOM,);
//    }
    //endregion

    dlPdlnbar = gsl_spline_eval (dlPdlnbarspline, h, dlPdlnbaracc);
}

// Depreciated methods TODO: remove when sure that read_cgs() and generate_data() work properly
void eos_table::read_cgs_old(string filename_in, double mu0, int print) {
    type    = "CGS tabulated EoS";
    filename = filename_in;

    //region Open file: ifstream fich
    ifstream fich(filename.c_str()) ;
    if (!fich) {
        cerr << "eos_table: " << "Can not open " << filename.c_str() <<"!"<< endl ;
        cerr << "Aborting..." << endl ;
        abort() ;
    }
    //endregion

    //region Process Header
    /* Header format example:
        # 2017.02.13 - 16:14
        # Polytropic EoS: (Kappa=0.05|Gamma=2)
        # M. J. Steil
        # Range: [9.992007E-15,1.718282E+00]fm**-3
        #
        80 <-- Number of lines
        #
        #        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
        #
    */
    fich.ignore(1000, '\n') ; // Jump over first line
    fich.ignore(1);
    safeGetline2(fich, info);// Get info (header line 2)

    for (int i=0; i<3; i++) {fich.ignore(1000, '\n');}; // Skip to line 6

    fich >> table_length ; fich.ignore(1000, '\n') ;   // Get table_length

    for (int i=0; i<3; i++) {fich.ignore(1000, '\n');}; // Skip column headers
    //endregion

    mB  = LmB_km;  // Mean baryon rest mass; standard is mBL, not used by this class directly

    //region Read in data
    vector<double> h_in(table_length);
    vector<double> lP_in(table_length);
    vector<double> dlP_in(table_length);

    vector<double> lnbar_in(table_length);
    vector<double> lrho_in(table_length);

    int no ;
    double nbar_fm3, rho_cgs, P_cgs;
    double nbar_tmp, rho_tmp, P_tmp, h_tmp;

    if(print)printf("-: i | lP_data[i], dlP_data[i], lrho_data[i], lnbar_data[i] | h_data[i], h_data[i]-h_data[i-1] \n");
    for (int i=0; i<table_length; i++) {
        fich >> no ;
        fich >> nbar_fm3 ;
        fich >> rho_cgs ;
        fich >> P_cgs ;
        fich.ignore(1000,'\n');

        // Interpolation is done for P and rho in CGS and for nbar in fm**-3 to be consistent with LORENE eos_tabul
        P_tmp = P_cgs/c0c0_cgs;   // [P_tmp] = gcm**-3
        nbar_tmp = nbar_fm3;    // [nbar_tmp] = fm**-3
        rho_tmp = rho_cgs;      // [rho_tmp] = gcm**-3

        // mu0/h0-Setup
        if(i==0) {
            if (mu0 == 0) {
                mu0scale = (P_tmp + rho_tmp) / (nbar_tmp);
            } else {
                mu0scale = mu0;
            }
            nbarScale = 1 / mu0scale * LrhoNuc_cgs * cfm3km3;
        }

        h_tmp = log((P_tmp+rho_tmp)/(nbar_tmp*mu0scale));

        h_in[i] = h_tmp;
        lP_in[i] = log10(P_tmp/LrhoNuc_cgs);
        dlP_in[i] = (rho_tmp+P_tmp)/P_tmp*h_tmp;

        lnbar_in[i] = log10(nbar_tmp);
        lrho_in[i] = log10(rho_tmp);

        if(i>0){
            if(print)printf("-: %d | %.4E, %.4E, %.4E, %.4E | %.4E, %d \n",i,lP_in[i],dlP_in[i],lrho_in[i],lnbar_in[i],h_in[i],(h_in[i]>h_in[i-1]));
        }
    }
    fich.close();

    h0 = h_in[0];
    P0 = exp10(lP_in[0])*LrhoNuc_km;
    rho0 = exp10(lrho_in[0])*cgcm3km2;
    nbar0 = exp10(lnbar_in[0])*cfm3km3;
    hmax = h_in[table_length-1];

    h_data.assign(h_in.begin(),h_in.begin()+table_length);
    lP_data.assign(lP_in.begin(),lP_in.begin()+table_length);
    dlP_data.assign(dlP_in.begin(),dlP_in.begin()+table_length);

    lnbar_data.assign(lnbar_in.begin(),lnbar_in.begin()+table_length);
    lrho_data.assign(lrho_in.begin(),lrho_in.begin()+table_length);
    //endregion

    //region Compute dlPdlnbar_data using cubic polynomials and setup dlPdlnbarspline
    vector<double> dlPdlnbar_in(table_length);
    double p0, p1, p2, n0, n1, n2, dpdnb;

    // special case: i=0 - Surface
    p0 = lP_in[0];
    p1 = lP_in[1];
    p2 = lP_in[2];

    n0 = lnbar_in[0];
    n1 = lnbar_in[1];
    n2 = lnbar_in[2];

    dpdnb = p0*(2*n0-n1-n2)/(n0-n1)/(n0-n2) +
            p1*(n0-n2)/(n1-n0)/(n1-n2) +
            p2*(n0-n1)/(n2-n0)/(n2-n1) ;

    dlPdlnbar_in[0]= dpdnb;

    // Main loop 1<=i<table_length-1
    for(int i=1;i<table_length-1;i++) {
        p0 = lP_in[i-1];
        p1 = lP_in[i];
        p2 = lP_in[i+1];

        n0 = lnbar_in[i-1];
        n1 = lnbar_in[i];
        n2 = lnbar_in[i+1];

        dpdnb = p0*(n1-n2)/(n0-n1)/(n0-n2) +
                p1*(2*n1-n0-n2)/(n1-n0)/(n1-n2) +
                p2*(n1-n0)/(n2-n0)/(n2-n1) ;

        dlPdlnbar_in[i]= dpdnb;
    }

    // special case: i=table_length-1 - Nucleus
    p0 = lP_in[table_length-3];
    p1 = lP_in[table_length-2];
    p2 = lP_in[table_length-1];

    n0 = lnbar_in[table_length-3];
    n1 = lnbar_in[table_length-2];
    n2 = lnbar_in[table_length-1];

    dpdnb = p0*(n2-n1)/(n0-n1)/(n0-n2) +
            p1*(n2-n0)/(n1-n0)/(n1-n2) +
            p2*(2*n2-n0-n1)/(n2-n0)/(n2-n1);

    dlPdlnbar_in[table_length-1]= dpdnb;

    // Assign data
    dlPdlnbar_data.assign(dlPdlnbar_in.begin(),dlPdlnbar_in.begin()+table_length);
    //endregion

    stringstream name_tmp;
    name_tmp << "Tabulated EoS (CGS):" << info.c_str() <<" - (from: " << filename.c_str() <<", "<< table_length <<" points)";
    name    = name_tmp.str();
}

