#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************/
/*Raw output is loged quantity.*/
/*Exponentiating must be done to be used outside this code*/
/************************/

#define AEQ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define PI 3.141592
#define Gammad -3
#define DNeff 0.6384 /*Delta N_eff at a_tr*/
#define ainit 1e-10
#define atr 5.10479e-5

/*Conversion factors*/
#define KG_to_GeV 5.62e26
#define M_to_GeV 5.07e15
#define S_to_GeV 1.519e24
#define M_to_Mpc 3.09e22
#define Gyr_to_Mpc 306.6
#define C 299792458
#define GeV_to_Mpc 1.563e28

#define Omega_ddm 2.78e-6
#define Omega_b 0.0443678
#define Omega_cdm 0.2523
#define Omega_Lambda 0.679
#define Omega_photon 5.4e-5
#define Omega_neutrino 1e-6

#define H0 1.44e-33 /*eV^-1*/
#define RHOC 3.67e-11 /*eV^4*/
#define G 6.70746e-57 /*eV^-2*/


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
/************************************************************************************* */

void Boltzmann(double x0, double ddm0, double CFT0, double h, int step, double Gamma0, double *DDM, double *CFT)
{
    FILE *output;
    output = fopen("output.dat", "w");
    double x = x0;
    double rho_ddm = log(ddm0);
    double rho_CFT = log(CFT0);

    DDM[0] = exp(rho_ddm);
    CFT[0] = exp(rho_CFT);

    fprintf(output, "# x\t\tlog(rho_ddm)\t\tlog(rho_CFT)\n");
    fprintf(output, "%e\t%e\t%e\n", x, rho_ddm, rho_CFT);

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

        fprintf(output, "%e\t%e\t%e\n", x, rho_ddm, rho_CFT);
    }
    fclose(output);
    printf("Result is saved in output.dat.\n");
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
    double h = -1e-4;
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

    Boltzmann(log(atr), rho_ddm_atr, rho_CFT_atr, h, N, Gamma0, DDM, CFT);

    free(DDM); free(CFT); free(a);

    return 0;
}