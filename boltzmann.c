#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AEQ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define HUBBLE 67.32 /*H_0 (km * s^-1 * Mpc^-1)*/
#define RHOC 8.505e-27 /*Today's critical density (kg * m^-3)*/
#define PI 3.141592
#define G 6.674e-11 /*Newton's constant (m^3 * kg^-1 * s^-2)*/

/*Conformal Hubble constant*/


int main()
{
    double h = -1e-3;
    double ainit = 1e-7;
    double atr = 1e-4; /*Scale factor at which decay rate is turned off*/
    double DNeff = 1.0; /*Delta N_eff at a_tr*/
    double Gamma0 = 1e-3, Gammad = 1.0; /*Gamma(a) ~ Gamma_0 * a^(Gamma_d) in Gyr^-1*/
    /*Today's Density parameters*/
    double Omega_ddm = 1e-5, Omega_b = 0.049389, Omega_cdm = 0.265028, Omega_Lambda = 0.679, Omega_photon = 1e-5; /*Flat universe*/
    double Omega_neutrino = 1.0 - (Omega_b + Omega_cdm + Omega_ddm + Omega_Lambda + Omega_photon);

    double Omega_CFT_atr = Omega_photon*pow(atr, -4)*(7/8)*pow(4/11, 4/3)*DNeff;
    double Omega_ddm_atr = Omega_ddm*pow(atr, -3);

    double *rho_ddm = (double *)malloc(sizeof(double)*(int)((log(atr) - log(ainit))/abs(h)));
    double *rho_CFT = (double *)malloc(sizeof(double)*(int)((log(atr) - log(ainit))/abs(h)));

    printf("%f\n",((log(atr) - log(ainit))/abs(h)));

    free(rho_ddm);
    free(rho_CFT);

    return 0;
}