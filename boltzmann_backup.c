#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/************************/
/*Raw output is loged quantity.*/
/*Exponentiating must be done to be used outside this code*/
/************************/

#define AEQ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define PI 3.141592
#define Gammad -3.5
#define DNeff 0.6384 /*Delta N_eff at a_tr*/
#define ainit 1e-10
#define atr 5.10479e-5

#define Omega_ddm 0.2523
#define Omega_b 0.0443678
#define Omega_cdm 0
#define Omega_Lambda 0.679
#define Omega_photon 5.4e-5
#define Omega_neutrino 1e-6

#define H0 1.44e-33 /*eV^-1*/
#define RHOC 3.67e-11 /*eV^4*/
#define G 6.70746e-57 /*eV^-2*/

int coeff[3] = {3.0/2.0, -2.0, 1.0/2.0};


/*Conformal Hubble constant in eV^-1*/
double H(double x, double rho_ddm, double rho_CFT)
{
    double rhotot;
    rhotot = rho_CFT + rho_ddm + RHOC*((Omega_b + Omega_cdm)*exp(-3*x) + (Omega_photon + Omega_neutrino)*exp(-4*x) + Omega_Lambda*exp(0.039*x));
    if (rhotot <= 0){
        printf("Negative energy at a = %e!\n", exp(x));
        exit(1);
    }
    double log_H2 = 2*x + log(rhotot) + log(8*PI*G/3.0);
    return exp(0.5*log_H2);
}

/*Evolution equation of log(rho)*/
double DDMEQ(double x, double rho_ddm, double rho_CFT, double Gamma0)
{
    return -3.0 - (Gamma0/H(x, exp(rho_ddm), exp(rho_CFT))) * exp((Gammad + 1)*x);
}

double CFTEQ(double x, double rho_ddm, double rho_CFT, double Gamma0)
{
    return -4.0 + (Gamma0/H(x, exp(rho_ddm), exp(rho_CFT)) * exp((Gammad + 1)*x))*exp(rho_ddm - rho_CFT);
}

void deriv(double x, double *rho, double *drho, double Gamma0)
{
    double h = H(x, rho[0], rho[1]);
    drho[0] = -3*rho[0] - (Gamma0/h)*exp((Gammad+1)*x)*rho[0];
    drho[1] = -4*rho[1] + (Gamma0/h)*exp((Gammad+1)*x)*rho[0];
    return;
}

void bdf2(double x_next, double h, double y_prev[], double y_prev2[], double y_next[])
{
    double y_guess[2];
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < 10; ++iter){
        double F[2];
        double df[2];
        deriv(x_next, y_guess, df);

        for (int i = 0; i < 2; ++i){
            F[i] = coeff[0]*y_guess[i] + coeff[1]*y_prev[i] + coeff[2]*y_prev2[i] - h*df[i];
            
        }
    }
}
/************************************************************************************* */

void Boltzmann(double x0, double ddm0, double CFT0, double h, int step, double Gamma0, double *DDM, double *CFT, char *s)
{
    FILE *output;
    output = fopen(s, "w");
    double x = x0;
    double rho_ddm = log(ddm0);
    double rho_CFT = log(CFT0);

    DDM[0] = exp(rho_ddm);
    CFT[0] = exp(rho_CFT);

    fprintf(output, "# x\t\tlog(rho_ddm)\t\tlog(rho_CFT)\n");
    fprintf(output, "%e\t%e\t%e\t%e\t%e\n", x, rho_ddm, rho_CFT, H(x, exp(rho_ddm), exp(rho_CFT)), Gamma0*exp(Gammad*x));

    for (int i = 0; i < step; i++){
        double k1_1 = DDMEQ(x, rho_ddm, rho_CFT, Gamma0);
        double k1_2 = CFTEQ(x, rho_ddm, rho_CFT, Gamma0);

        double k2_1 = DDMEQ(x+h/2.0, rho_ddm+h*k1_1/2.0, rho_CFT+h*k1_2/2.0, Gamma0);
        double k2_2 = CFTEQ(x+h/2.0, rho_ddm+h*k1_1/2.0, rho_CFT+h*k1_2/2.0, Gamma0);

        double k3_1 = DDMEQ(x+h/2.0, rho_ddm+h*k2_1/2.0, rho_CFT+h*k2_2/2.0, Gamma0);
        double k3_2 = CFTEQ(x+h/2.0, rho_ddm+h*k2_1/2.0, rho_CFT+h*k2_2/2.0, Gamma0);

        double k4_1 = DDMEQ(x+h, rho_ddm+h*k3_1, rho_CFT+h*k3_2, Gamma0);
        double k4_2 = CFTEQ(x+h, rho_ddm+h*k3_1, rho_CFT+h*k3_2, Gamma0);

        rho_ddm += h*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6.0;
        rho_CFT += h*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6.0;
        x += h;

        DDM[i+1] = exp(rho_ddm);
        CFT[i+1] = exp(rho_CFT);

        fprintf(output, "%e\t%e\t%e\t%e\t%e\n", x, rho_ddm, rho_CFT, H(x, exp(rho_ddm), exp(rho_CFT)), Gamma0*exp(Gammad*x));
    }
    fclose(output);
    printf("Result is saved in %s.\n", s);
}

/*Interpolation after atr*/
/*rho(a) = rho(atr)*a^-3 (a > atr)*/
double interpolation(double *rho, double a)
{
    if (a >= atr){
        return rho[0]*pow(atr/a, 3);
    }
    else{
        printf("Target a must be larger than atr/\n");
        exit(1);
    }
}

int main()
{
    /*General parameters*/
    double h = -1e-5;
    double Gamma0 = 1e3; /*Gamma(a) ~ Gamma_0 * a^(Gamma_d) in Gyr^-1*/
    int N = (int)((log(atr) - log(ainit))/fabs(h)); /*The number of steps*/

    /*Convert to eV units*/
    Gamma0 *= 2.085e-50; /*eV^-1*/

    double rho_CFT_atr = RHOC*Omega_photon*pow(atr, -4.0)*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0)*DNeff;
    double rho_ddm_atr = RHOC*Omega_ddm*pow(atr, -3.0);
    
    /*Original data (not log)*/
    double *DDM = (double *)malloc(sizeof(double)*(N+1));
    double *CFT = (double *)malloc(sizeof(double)*(N+1));
    double *a = (double *)malloc(sizeof(double)*(N+1));

    for (int i = 0; i < N; i++){
        a[i] = log(atr) + i*h;
        a[i] = exp(a[i]);
    }

    Boltzmann(log(atr), rho_ddm_atr, rho_CFT_atr, h, N, Gamma0, DDM, CFT, "output_backward.dat");
    Boltzmann(log(ainit), DDM[N], CFT[N], -1*h, N, Gamma0, DDM, CFT, "output_forward.dat");
    FILE *params;
    params = fopen("params.dat", "w");
    fprintf(params, "%s\t%e\n", "aeq", AEQ);
    fprintf(params, "%s\t%e\n", "Gammad", Gammad);
    fprintf(params, "%s\t%e\n", "Gamma0", Gamma0);
    fprintf(params, "%s\t%e\n", "DNeff", DNeff);
    fprintf(params, "%s\t%e\n", "atr", atr);
    fclose(params);

    free(DDM); free(CFT); free(a);

    return 0;
}