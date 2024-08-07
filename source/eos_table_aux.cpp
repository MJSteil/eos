//
// Created by M. J. Steil on 2017.07.04.
//

#include "../eos_table.hpp"
#include "../../matplotlib-cpp/matplotlibcpp.h"

// Auxillary
void eos_table::plot_interpolation(int plot) {

    int ij=0;
    vector<double> grid_h, grid_P, grid_rho, grid_nbar,grid_drho,grid_dP,grid_dnbar;
    vector<double> lhi,hi, Pi,rhoi,nbari,drhoi,dPi,dnbari,csi;

    for(int i=0;i<table_length-1;i++){
        grid_h.push_back(exp10(lh_data[i]));

        grid_P.push_back(exp10(lP_data[i])*LrhoNuc_km/cMeVfm3km2);
        grid_rho.push_back(exp10(lrho_data[i])*cgcm3km2/cMeVfm3km2);
        grid_nbar.push_back(exp10(lnbar_data[i]));

        grid_dP.push_back(dlPdlh(exp10(lh_data[i])));
        grid_drho.push_back(dlrhodlh(exp10(lh_data[i])));
        grid_dnbar.push_back(dlnbardlh(exp10(lh_data[i])));

        for(int j=0;j<9;j++){
            if(i<table_length-1){
                double h = lh_data[i]+(lh_data[i+1]-lh_data[i])/9*j;
                hi.push_back(exp10(h));

                Pi.push_back(P(exp10(h))/cMeVfm3km2);
                rhoi.push_back(rho(exp10(h))/cMeVfm3km2);
                nbari.push_back(nbar(exp10(h))/cfm3km3);

                dPi.push_back(dlPdlh(exp10(h)));
                drhoi.push_back(dlrhodlh(exp10(h)));
                dnbari.push_back(dlnbardlh(exp10(h)));

                csi.push_back(cs(h,0));
                ij++;

            }

        }
    }
    switch ( plot )
    {
        case 1:
            matplotlibcpp::loglog(hi,Pi);
            matplotlibcpp::loglog(grid_h,grid_P,"b.");
            matplotlibcpp::title("P over h - 1");
            matplotlibcpp::xlabel("h");
            matplotlibcpp::ylabel("P (MeV/fm^-3)");
            matplotlibcpp::show();
            break;
        case 2:
            matplotlibcpp::loglog(hi,rhoi);
            matplotlibcpp::loglog(grid_h,grid_rho,"b.");
            matplotlibcpp::title("rho over h - 2");
            matplotlibcpp::xlabel("h");
            matplotlibcpp::ylabel("rho (MeV/fm^-3)");
            matplotlibcpp::show();
            break;
        case 3:
            matplotlibcpp::loglog(hi,nbari);
            matplotlibcpp::loglog(grid_h,grid_nbar,"b.");
            matplotlibcpp::title("nbar over h - 3");
            matplotlibcpp::xlabel("h");
            matplotlibcpp::ylabel("n (fm^-3)");
            matplotlibcpp::show();
            break;
        case 4:
            matplotlibcpp::semilogx(hi,csi);
            matplotlibcpp::title("cs over h - 3");
            matplotlibcpp::xlabel("h");
            matplotlibcpp::ylabel("cs");
            matplotlibcpp::show();
            break;
        default:
            matplotlibcpp::loglog(rhoi,Pi);
            matplotlibcpp::loglog(grid_rho,grid_P,"b.");
            matplotlibcpp::title("P over rho - 1");
            matplotlibcpp::xlabel("rho (MeV/fm^-3)");
            matplotlibcpp::ylabel("P (MeV/fm^-3)");
            matplotlibcpp::show();
    }
}