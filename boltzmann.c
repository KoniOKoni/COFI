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
#define Gammad -5.0
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

#define TOL 1e-10
#define MAX_ITER 100

double coeff[3] = {3.0/2.0, -2.0, 1.0/2.0};

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

/*Jacobian*/
void jacobian(double x, double *rho, double Gamma0 ,double jac[][2], double a, double h)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            jac[i][j] = 0;
        }
    }
    double Gamma = Gamma0*exp((Gammad + 1)*x);
    double Hub = H(x, exp(rho[0]), exp(rho[1]));
    double k = Gamma/Hub;

    jac[0][0] = (a - exp(rho[1] - rho[0])*h*k)/(pow(a,2) - a*exp(rho[1]-rho[0])*h*k);
    jac[0][1] = 0.0;
    jac[1][0] = -1*exp(rho[1]-rho[0])*h*k/(pow(a,2) - a*exp(rho[1]-rho[0])*h*k);
    jac[1][1] = a/(pow(a,2) - a*exp(rho[1]-rho[0])*h*k);
}

/*Evolution equation of rho*/
void deriv(double x, double *rho, double *drho, double Gamma0)
{
    double h = H(x, exp(rho[0]), exp(rho[1]));
    drho[0] = -3 - (Gamma0/h)*exp((Gammad+1)*x);
    drho[1] = -4 + (Gamma0/h)*exp((Gammad+1)*x)*exp(rho[1] - rho[0]);
}

void bdf2(double x_next, double h, double Gamma0 ,double y_prev[], double y_prev2[], double y_next[], double jac[][2])
{
    double y_guess[2];
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < MAX_ITER; iter++){
        double F[2];
        double df[2];
        double delta[2] = {0.0,0.0};
        deriv(x_next, y_guess, df, Gamma0);
        jacobian(x_next, y_guess, Gamma0, jac, coeff[0], h);

        for (int i = 0; i < 2; i++){
            F[i] = coeff[0]*y_guess[i] + coeff[1]*y_prev[i] + coeff[2]*y_prev2[i] - h*df[i];
            for (int j = 0; j < 2; j++){
                delta[i] += -F[j]*jac[i][j];
            }
            y_guess[i] += delta[i];
        }

        double max_delta = 0.0;
        for (int i = 0; i < 2; i++){
            if (fabs(F[i]) > max_delta) max_delta = fabs(F[i]); 
        }
        if (max_delta < TOL) break;
    }
    memcpy(y_next, y_guess, sizeof(double) * 2);
}
/************************************************************************************* */

/*Interpolation*/
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
    double h = 1e-5;
    double Gamma0 = 1e3; /*Gamma(a) ~ Gamma_0 * a^(Gamma_d) in Gyr^-1*/
    int N = (int)((log(atr) - log(ainit))/fabs(h)); /*The number of steps*/

    /*Convert to eV units*/
    Gamma0 *= 2.085e-50; /*eV^-1*/

    /*double rho_CFT_atr = RHOC*Omega_photon*pow(atr, -4.0)*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0)*DNeff;
    double rho_ddm_atr = RHOC*Omega_ddm*pow(atr, -3.0);*/
    double rho_CFT_atr = 1.73e24;
    double rho_ddm_atr = 1.58e60;
    
    /*Original data (not log)*/
    double jac[2][2] = {{0.0,0.0},{0.0,0.0}};

    double x = log(ainit);
    double rho_prev2[2] = {log(rho_ddm_atr) - 3*h, log(rho_CFT_atr) - 4*h};
    double rho_prev[2];
    double rho_next[2];

    FILE *output, *params;
    output = fopen("output.dat", "w");

    memcpy(rho_prev, rho_prev2, sizeof(double)*2);
    for (int iter = 0; iter < MAX_ITER; iter++){
        double df[2];
        double F[2];
        double delta[2] = {0.0, 0.0};
        deriv(x+h, rho_prev, df, Gamma0);
        jacobian(x+h, rho_prev, Gamma0, jac, 1, h);
        for (int i = 0; i < 2; i++){
            F[i] = rho_prev[i] - rho_prev2[i] - h*df[i];
            for (int j = 0; j < 2; j++){
                delta[i] += -F[j]*jac[i][j];
            }
            rho_prev[i] += delta[i];
        }
    }

    for (int i = 0; i < N; i++){
        x += h;
        bdf2(x, h, Gamma0, rho_prev, rho_prev2, rho_next, jac);
        fprintf(output, "%e\t%e\t%e\t%e\t%e\n", x, rho_next[0], rho_next[1], H(x, exp(rho_next[0]), exp(rho_next[1])), Gamma0*exp((Gammad + 1)*x));

        memcpy(rho_prev2, rho_prev, sizeof(double) * 2);
        memcpy(rho_prev, rho_next, sizeof(double) * 2);
    }

    params = fopen("params.dat", "w");
    fprintf(params, "%s\t%e\n", "aeq", AEQ);
    fprintf(params, "%s\t%e\n", "Gammad", Gammad);
    fprintf(params, "%s\t%e\n", "Gamma0", Gamma0);
    fprintf(params, "%s\t%e\n", "DNeff", DNeff);
    fprintf(params, "%s\t%e\n", "atr", atr);
    fprintf(params, "%s\t%e\n", "ainit", ainit);
    fclose(params); fclose(output);

    return 0;
}