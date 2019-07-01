#include <math.h>
#include <stdio.h>
#include <float.h>
#include "mtrand.h"
#define exponentiale M_E
#define pi M_PI
double max(double a, double b){
return a > b ? a : b;}
double min(double a, double b){
return a < b ? a : b;}
double root(double n,double x);

double root_0(double n,double x);

double root_1(double n,double x);

double cot(double x);

double cot_0(double x);

double arccot(double x);

double arccot_0(double x);

double coth(double x);

double coth_0(double x);

double csc(double x);

double csc_0(double x);

double arccsc(double x);

double arccsc_0(double x);

double csch(double x);

double csch_0(double x);

double sec(double x);

double sec_0(double x);

double arcsec(double x);

double arcsec_0(double x);

double sech(double x);

double sech_0(double x);

void res_function_(double *time_ptr, double *dynamicVars, double *yprime, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar);

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res);

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual);

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *pd, double *cj_ptr, double *constants, int *intpar);

void dres_dkloffPI3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd);

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar);

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar);

void integrate_stochatic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory);

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar);

double root(double n,double x){
return pow(x, 1.0/n);
}

double root_0(double n,double x){
return -(log(x)*pow(x, 1.0/n)*1.0/pow(n, 2.0));
}

double root_1(double n,double x){
return pow(x, 1.0/n - 1.0)/n;
}

double cot(double x){
return 1.0/tan(x);
}

double cot_0(double x){
return -(1.0/(pow(cos(x), 2.0)*pow(tan(x), 2.0)));
}

double arccot(double x){
return atan(1.0/x);
}

double arccot_0(double x){
return -(1.0/pow(x, 2.0)/(pow(1.0/x, 2.0) + 1.0));
}

double coth(double x){
return 1.0/tanh(x);
}

double coth_0(double x){
return -(1.0/(pow(cosh(x), 2.0)*pow(tanh(x), 2.0)));
}

double csc(double x){
return 1.0/sin(x);
}

double csc_0(double x){
return -(cos(x)/pow(sin(x), 2.0));
}

double arccsc(double x){
return asin(1.0/x);
}

double arccsc_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double csch(double x){
return 1.0/sinh(x);
}

double csch_0(double x){
return -(cosh(x)/pow(sinh(x), 2.0));
}

double sec(double x){
return 1.0/cos(x);
}

double sec_0(double x){
return sin(x)/pow(cos(x), 2.0);
}

double arcsec(double x){
return acos(1.0/x);
}

double arcsec_0(double x){
return -(1.0/pow(x, 2.0)/sqrt(1.0 - pow(1.0/x, 2.0)));
}

double sech(double x){
return 1.0/cosh(x);
}

double sech_0(double x){
return -(sinh(x)/pow(cosh(x), 2.0));
}

void res_function_(double *time_ptr, double *dynamicVars, double *yprime, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

residual[0] = -yprime[0];
residual[1] = -yprime[1];
residual[2] = -yprime[2];
residual[3] = -yprime[3];
residual[4] = kc*(1.0 - switch0)*circ*3.0 + kcoff*probeCLTA + kcoff*probeCLTAPI3P - ku*CLTAfree - kcon*Nprobe*CLTAfree - kconb*probePI3P*CLTAfree - yprime[4];
residual[5] = -yprime[5];
residual[6] = -yprime[6];
residual[7] = -(kcat*switch0*N3k*PI/(PI + KM)) - yprime[7];
residual[8] = kcat*switch0*N3k*PI/(PI + KM) + kloffPI3P*probePI3P + kloffPI3P*probeCLTAPI3P - kcat*switch0*N4k*PI3P/(PI3P + KM) - klon*Np*PI3P - klonb*probeCLTA*PI3P - yprime[8];
residual[9] = kcat*switch0*N5p*PI45P2/(PI45P2 + KM) - kcat*switch0*N3k4*PI4P/(PI4P + KM) - yprime[9];
residual[10] = -(kcat*switch0*N5p*PI45P2/(PI45P2 + KM)) - yprime[10];
residual[11] = kcat*switch0*N4k*PI3P/(PI3P + KM) + kcat*switch0*N3k4*PI4P/(PI4P + KM) - yprime[11];
residual[12] = kcon*Nprobe*CLTAfree + kloffPI3P*probeCLTAPI3P - ku*probeCLTA - klonb*probeCLTA*PI3P - kcoff*probeCLTA - yprime[12];
residual[13] = ku*probeCLTAPI3P + klon*Np*PI3P + kcoff*probeCLTAPI3P - kconb*probePI3P*CLTAfree - kloffPI3P*probePI3P - yprime[13];
residual[14] = klonb*probeCLTA*PI3P + kconb*probePI3P*CLTAfree - ku*probeCLTAPI3P - kloffPI3P*probeCLTAPI3P - kcoff*probeCLTAPI3P - yprime[14];
}

void alg_deriv_func_(double *alg_yp, double *dynamicVars, double *yp, double *time_ptr, double *constants, double *alg_derivs_res){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

}

void alg_res_func_(double *alg_vals, double *dynamicVars, double *time_ptr, double *constants, double *residual){
double time = *time_ptr;


double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

}

void dres_dc_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

pd[64] = kc*(1.0 - switch0)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*3.0 - ku - CLTAfree*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333) - kcon*Nprobe - kconb*probePI3P;
pd[79] = -(kc*circ*3.0) - CLTAfree*(circ + fleak*CLTA)*(1.0 - switch1)*(ku1 - ku0)/(CLTA + 0.01);
pd[94] = -(CLTAfree*(circ + fleak*CLTA)*(ku2 - ((1.0 - switch0)*ku0 + switch0*ku1))/(CLTA + 0.01));
pd[184] = kcoff + kc*(1.0 - switch0)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*3.0 - CLTAfree*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[199] = -(kconb*CLTAfree);
pd[214] = kcoff + kc*(1.0 - switch0)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*3.0 - CLTAfree*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[7] = -(kcat*switch0*PI/(PI + KM));
pd[82] = -(kcat*N3k*PI/(PI + KM));
pd[112] = -(kcat*switch0*N3k/(PI + KM) - kcat*switch0*N3k*PI/pow(PI + KM, 2.0));
pd[8] = kcat*switch0*PI/(PI + KM);
pd[23] = -(kcat*switch0*PI3P/(PI3P + KM));
pd[83] = kcat*N3k*PI/(PI + KM) - kcat*N4k*PI3P/(PI3P + KM);
pd[113] = kcat*switch0*N3k/(PI + KM) - kcat*switch0*N3k*PI/pow(PI + KM, 2.0);
pd[128] = kcat*switch0*N4k*PI3P/pow(PI3P + KM, 2.0) - kcat*switch0*N4k/(PI3P + KM) - klon*Np - klonb*probeCLTA;
pd[188] = -(klonb*PI3P);
pd[203] = kloffPI3P;
pd[218] = kloffPI3P;
pd[39] = kcat*switch0*PI45P2/(PI45P2 + KM);
pd[54] = -(kcat*switch0*PI4P/(PI4P + KM));
pd[84] = kcat*N5p*PI45P2/(PI45P2 + KM) - kcat*N3k4*PI4P/(PI4P + KM);
pd[144] = kcat*switch0*N3k4*PI4P/pow(PI4P + KM, 2.0) - kcat*switch0*N3k4/(PI4P + KM);
pd[159] = kcat*switch0*N5p/(PI45P2 + KM) - kcat*switch0*N5p*PI45P2/pow(PI45P2 + KM, 2.0);
pd[40] = -(kcat*switch0*PI45P2/(PI45P2 + KM));
pd[85] = -(kcat*N5p*PI45P2/(PI45P2 + KM));
pd[160] = -(kcat*switch0*N5p/(PI45P2 + KM) - kcat*switch0*N5p*PI45P2/pow(PI45P2 + KM, 2.0));
pd[26] = kcat*switch0*PI3P/(PI3P + KM);
pd[56] = kcat*switch0*PI4P/(PI4P + KM);
pd[86] = kcat*N4k*PI3P/(PI3P + KM) + kcat*N3k4*PI4P/(PI4P + KM);
pd[131] = kcat*switch0*N4k/(PI3P + KM) - kcat*switch0*N4k*PI3P/pow(PI3P + KM, 2.0);
pd[146] = kcat*switch0*N3k4/(PI4P + KM) - kcat*switch0*N3k4*PI4P/pow(PI4P + KM, 2.0);
pd[72] = kcon*Nprobe - probeCLTA*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[87] = -(probeCLTA*(circ + fleak*CLTA)*(1.0 - switch1)*(ku1 - ku0)/(CLTA + 0.01));
pd[102] = -(probeCLTA*(circ + fleak*CLTA)*(ku2 - ((1.0 - switch0)*ku0 + switch0*ku1))/(CLTA + 0.01));
pd[132] = -(klonb*probeCLTA);
pd[192] = -ku - probeCLTA*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333) - klonb*PI3P - kcoff;
pd[222] = kloffPI3P - probeCLTA*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[73] = probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333) - kconb*probePI3P;
pd[88] = probeCLTAPI3P*(circ + fleak*CLTA)*(1.0 - switch1)*(ku1 - ku0)/(CLTA + 0.01);
pd[103] = probeCLTAPI3P*(circ + fleak*CLTA)*(ku2 - ((1.0 - switch0)*ku0 + switch0*ku1))/(CLTA + 0.01);
pd[133] = klon*Np;
pd[193] = probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[208] = -(kconb*CLTAfree) - kloffPI3P;
pd[223] = ku + kcoff + probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[74] = kconb*probePI3P - probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[89] = -(probeCLTAPI3P*(circ + fleak*CLTA)*(1.0 - switch1)*(ku1 - ku0)/(CLTA + 0.01));
pd[104] = -(probeCLTAPI3P*(circ + fleak*CLTA)*(ku2 - ((1.0 - switch0)*ku0 + switch0*ku1))/(CLTA + 0.01));
pd[134] = klonb*probeCLTA;
pd[194] = klonb*PI3P - probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333);
pd[209] = kconb*CLTAfree;
pd[224] = -ku - probeCLTAPI3P*(((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0))*(1.0/CLTAmax - CLTA*2.0/pow(CLTAmax, 2.0))*0.166666666667/(sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))*sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)))/(CLTA + 0.01) + (((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*fleak/(CLTA + 0.01) - ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/pow(CLTA + 0.01, 2.0))*0.333333333333) - kloffPI3P - kcoff;
}

void dres_dcdot_function_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

pd[0] = -1;
pd[16] = -1;
pd[32] = -1;
pd[48] = -1;
pd[64] = -1;
pd[80] = -1;
pd[96] = -1;
pd[112] = -1;
pd[128] = -1;
pd[144] = -1;
pd[160] = -1;
pd[176] = -1;
pd[192] = -1;
pd[208] = -1;
pd[224] = -1;
}

void ddaskr_jac_(double *time_ptr, double *dynamicVars, double *yprime, double *pd, double *cj_ptr, double *constants, int *intpar){
double cj = *cj_ptr;

dres_dc_function_(time_ptr, dynamicVars, yprime, constants, pd);

double local_dres_dcdot[15*15] = {0};
dres_dcdot_function_(time_ptr, dynamicVars, yprime, constants, local_dres_dcdot);

int ii;
for(ii=0; ii < 225; ii++){
  pd[ii] += cj*local_dres_dcdot[ii];}
}

void dres_dkloffPI3P_(double *time_ptr, double *dynamicVars, double *yprime, double *constants, double *pd){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N4k = dynamicVars[1];
double N5p = dynamicVars[2];
double N3k4 = dynamicVars[3];
double CLTAfree = dynamicVars[4];
double switch0 = dynamicVars[5];
double switch1 = dynamicVars[6];
double PI = dynamicVars[7];
double PI3P = dynamicVars[8];
double PI4P = dynamicVars[9];
double PI45P2 = dynamicVars[10];
double PI34P2 = dynamicVars[11];
double probeCLTA = dynamicVars[12];
double probePI3P = dynamicVars[13];
double probeCLTAPI3P = dynamicVars[14];

double kc = ku0 + knet;
double kcon = kcon0/(CLTAmax*3.0);
double kconb = kcon/v_close;
double klon = klon_max/N45;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double klonb = klon/v_close;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;

pd[8] = probePI3P + probeCLTAPI3P;
pd[12] = probeCLTAPI3P;
pd[13] = -probePI3P;
pd[14] = -probeCLTAPI3P;
}

void sens_rhs_(double *time_ptr, double *sens_y, double *sens_yp, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){

res_function_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);

int p_index = (int)constants[26];
double constants_only[26];
int jj;
for (jj = 0; jj < 26; jj++){
constants_only[jj] = constants[jj];}
double *dc_dp = &sens_y[15];
double *dcdot_dp = &sens_yp[15];
double *local_dres_dp = &sens_res[15];
int ii;
for(ii = 0; ii < 15; ii++){
local_dres_dp[ii] = 0;}
switch(p_index)
{
case 0 : dres_dkloffPI3P_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);
break;
}
double local_dres_dc[225] = {0};
dres_dc_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dc);
int row, col;
for(row = 0; row < 15; row++){
for(col = 0; col < 15; col++){
sens_res[row+15] += local_dres_dc[row + col*15]*dc_dp[col];}}
double local_dres_dcdot[225] = {0};
dres_dcdot_function_(time_ptr, sens_y, sens_yp, constants, local_dres_dcdot);
for(row = 0; row < 15; row++){
for(col = 0; col < 15; col++){
sens_res[row+15] += local_dres_dcdot[row + col*15]*dcdot_dp[col];}}
}

void res_function_logdv_(double *time_ptr, double *log_dv, double *log_yp, double *cj_ptr, double *residual, int *ires_ptr, double *constants, int *ipar){
double dynamicVars[15];
double yprime[15];
int ii;
for(ii = 0; ii < 15; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
res_function_(time_ptr, dynamicVars, yprime, cj_ptr, residual, ires_ptr, constants, ipar);
}

void root_func_logdv_(int *neq_ptr, double *time_ptr, double *log_dv, double *log_yp, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double dynamicVars[15];
double yprime[15];
int ii;
for(ii = 0; ii < 15; ii++){
dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);
yprime[ii] = log_yp[ii] * dynamicVars[ii];}
root_func_(neq_ptr, time_ptr, dynamicVars, yprime, nrt_ptr, root_devs, constants, ipar);
}

void sens_rhs_logdv_(double *time_ptr, double *sens_y_log, double *sens_yp_log, double *cj_ptr, double *sens_res, int *ires_ptr, double *constants, int *ipar){
double sens_y[30];
double sens_yp[30];
int ii;
for(ii = 0; ii < 15; ii++){
sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);
sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}
for(ii = 15; ii < 30; ii++){
sens_y[ii] = sens_y_log[ii];
sens_yp[ii] = sens_yp_log[ii];}
sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, ires_ptr, constants, ipar);
}

void integrate_stochastic_tidbit_(unsigned long* seed_ptr, int* reseed, double* time_ptr, int* dv, double* cv, double* rmsd_ptr, double* stop_time_ptr, double* trajectory) {
return;}

void root_func_(int *neq_ptr, double *time_ptr, double *dynamicVars, double *yprime, int *nrt_ptr, double *root_devs, double *constants, int *ipar){
double time = *time_ptr;

double pit = constants[0];
double NC = constants[1];
double N4 = constants[2];
double N45 = constants[3];
double NPI = constants[4];
double Nfree = constants[5];
double v_close = constants[6];
double CLTAmax = constants[7];
double kcat = constants[8];
double KM = constants[9];
double N3k0 = constants[10];
double N4k0 = constants[11];
double N5p0 = constants[12];
double N3k40 = constants[13];
double knet = constants[14];
double ku0 = constants[15];
double ku1 = constants[16];
double ku2 = constants[17];
double fleak = constants[18];
double tsciss = constants[19];
double delay = constants[20];
double kcon0 = constants[21];
double kcoff = constants[22];
double klon_max = constants[23];
double Nprobe = constants[24];
double kloffPI3P = constants[25];

double N3k = dynamicVars[0];
double N3k_deriv_wrt_time = yprime[0];
double N4k = dynamicVars[1];
double N4k_deriv_wrt_time = yprime[1];
double N5p = dynamicVars[2];
double N5p_deriv_wrt_time = yprime[2];
double N3k4 = dynamicVars[3];
double N3k4_deriv_wrt_time = yprime[3];
double CLTAfree = dynamicVars[4];
double CLTAfree_deriv_wrt_time = yprime[4];
double switch0 = dynamicVars[5];
double switch0_deriv_wrt_time = yprime[5];
double switch1 = dynamicVars[6];
double switch1_deriv_wrt_time = yprime[6];
double PI = dynamicVars[7];
double PI_deriv_wrt_time = yprime[7];
double PI3P = dynamicVars[8];
double PI3P_deriv_wrt_time = yprime[8];
double PI4P = dynamicVars[9];
double PI4P_deriv_wrt_time = yprime[9];
double PI45P2 = dynamicVars[10];
double PI45P2_deriv_wrt_time = yprime[10];
double PI34P2 = dynamicVars[11];
double PI34P2_deriv_wrt_time = yprime[11];
double probeCLTA = dynamicVars[12];
double probeCLTA_deriv_wrt_time = yprime[12];
double probePI3P = dynamicVars[13];
double probePI3P_deriv_wrt_time = yprime[13];
double probeCLTAPI3P = dynamicVars[14];
double probeCLTAPI3P_deriv_wrt_time = yprime[14];

double kc = ku0 + knet;
double kc_deriv_wrt_time = 0.0;
double kcon = kcon0/(CLTAmax*3.0);
double kcon_deriv_wrt_time = 0.0;
double kconb = kcon/v_close;
double kconb_deriv_wrt_time = 0.0;
double klon = klon_max/N45;
double klon_deriv_wrt_time = 0.0;
double CLTA = (CLTAfree + probeCLTA + probeCLTAPI3P)/3.0;
double CLTA_deriv_wrt_time = CLTAfree_deriv_wrt_time*0.333333333333 + probeCLTA_deriv_wrt_time*0.333333333333 + probeCLTAPI3P_deriv_wrt_time*0.333333333333;
double circ = sqrt(sqrt(pow(CLTA/CLTAmax - pow(CLTA/CLTAmax, 2.0), 2.0)));
double circ_deriv_wrt_time = 0.0;
double ku = ((1.0 - switch1)*((1.0 - switch0)*ku0 + switch0*ku1) + switch1*ku2)*(circ + fleak*CLTA)/(CLTA + 0.01);
double ku_deriv_wrt_time = (circ + fleak*CLTA)*(ku2 - ((1.0 - switch0)*ku0 + switch0*ku1))*switch1_deriv_wrt_time/(CLTA + 0.01) + (circ + fleak*CLTA)*(1.0 - switch1)*(ku1 - ku0)*switch0_deriv_wrt_time/(CLTA + 0.01);
double klonb = klon/v_close;
double klonb_deriv_wrt_time = 0.0;
double probe = probeCLTA + probePI3P + probeCLTAPI3P;
double probe_deriv_wrt_time = probeCLTA_deriv_wrt_time + probePI3P_deriv_wrt_time + probeCLTAPI3P_deriv_wrt_time;

root_devs[0] = (time > tsciss) - 0.5;
root_devs[1] = (time > (tsciss + delay)) - 0.5;
}
