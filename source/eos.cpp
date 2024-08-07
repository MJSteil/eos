//
// Created by M. J. Steil on 2017.02.08.
//

#include "../eos.hpp"
eos::eos() {

}

void eos::compare(eos &eos_in, double h) {
    printf("|-> Compare %s vs %s @ h=%.4E \n",name.c_str(),eos_in.name.c_str(),h);
    printf("|-> %s vs %s\n|\n",type.c_str(),eos_in.type.c_str());

    vector<double> points = {
            eos_in.nbar(h),
            eos_in.P(h),
            eos_in.rho(h),
            eos_in.dlnbardlh(h),
            eos_in.dlPdlh(h),
            eos_in.dlrhodlh(h)
    };

    vector<double> errors = {
            fabs(1-nbar(h)/eos_in.nbar(h)),
            fabs(1-P(h)/eos_in.P(h)),
            fabs(1-rho(h)/eos_in.rho(h)),
            fabs(1-dlnbardlh(h)/eos_in.dlnbardlh(h)),
            fabs(1-dlPdlh(h)/eos_in.dlPdlh(h)),
            fabs(1-dlrhodlh(h)/eos_in.dlrhodlh(h))
    };

    double mean, var;
    mean = gsl_stats_mean(errors.data(),1,6);
    var = gsl_stats_variance(errors.data(),1,6);

    printf("|=> %.4E = |1-nbar/nbar_in| (nbar_in=%.4E km**-3) \n",errors[0],points[0]);
    printf("|=> %.4E = |1-P/P_in| (P_in=%.4E km**-2) \n",errors[1],points[1]);
    printf("|=> %.4E = |1-rho/rho_in| (rho_in=%.4E km**-2) \n|\n",errors[2],points[2]);

    printf("|=> %.4E = |1-dlnbardlh/dlnbardlh_in| (dlnbardlh_in=%.4E km**-3) \n",errors[3],points[3]);
    printf("|=> %.4E = |1-dlPdlh/dlPdlh_in| (dlPdlh_in=%.4E) \n",errors[4],points[4]);
    printf("|=> %.4E = |1-dlrhodlh/dlrhodlh_in| (dlrhodlh_in=%.4E) \n|\n",errors[5],points[5]);

    printf("|=> %.4E = mean_error (var=%.4E) \n",mean,var);
}

void eos::info() {
    printf("%s\n",name.c_str());
    printf("|-> h0=%.8E hmax=%.8E\n",h0,hmax);
    printf("|--> P0=%.8E MeVfm**-3 Pmax=%.8E MeVfm**-3\n",P(h0)/cMeVfm3km2,P(hmax)/cMeVfm3km2);
    printf("|--> rho0=%.8E MeVfm**-3 rhomax=%.8E MeVfm**-3\n",rho(h0)/cMeVfm3km2,rho(hmax)/cMeVfm3km2);
    printf("|--> nbar0=%.8E fm**-3 nbarmax=%.8E  fm**-3 \n",nbar(h0)/cfm3km3,nbar(hmax)/cfm3km3);
    printf("|--> csmax=%.8E \n",cs(hmax));
}

void eos::nbar_sat( double h) {

    cout<<  nbar(h)/cfm3km3 << endl;
    cout<<(rho(h)/nbar(h)/cMeVkm-(mp+mn)/2.) << endl;
    cout<<  9*(dlPdlh(h)-dlnbardlh(h))/dlnbardlh(h)*P(h)/nbar(h)/cMeVkm << endl;

    cout<<  -1+ dlrhodlh(h)/dlnbardlh(h) << endl;
}

double eos::cs(double h, double hthM) {
    if(h>h0*hthM){
        return sqrt(fabs(dlPdlh(h)*P(h)/(dlrhodlh(h)*rho(h))));
    }else{
        return 0;
    }
}
