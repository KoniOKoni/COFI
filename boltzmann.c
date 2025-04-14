#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define _SUCCESS_ 0
#define class_alloc(pointer, size, error_message)                    \
  if (((pointer) = malloc(size)) == NULL) {                          \
    sprintf(error_message, "Could not allocate memory");             \
    return 1;                                                        \
  }

/************************/
/*Raw output is loged quantity.*/
/*Exponentiating must be done to be used outside this code*/
/************************/

#define _Mpc_over_m_ 3.085677581282e22  /**< conversion factor from meters to megaparsecs */
#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears
                                       (c=1 units, Julian years of 365.25 days) */
#define _c_ 2.99792458e8            /**< c in m/s */
#define _G_ 6.67428e-11             /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */

#define _AEQ_ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define _PI_ 3.141592
#define _AINIT_ 1e-10

#define _TOL_ 1e-20
#define _MAX_ITER_ 100
#define _H_STEP_ 0.1
struct background {
    double *a_table;
    double *rho_chi_table;
    double *rho_cft_table;
    double *Gamma_table;
    double *H_table;
    double jac[2][2];
    int a_size;
    char error_message[1024];

    /*Today's value*/
    double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */
    double Omega0_g; /**< \f$ \Omega_{0 \gamma} \f$: photons */
    double Omega0_b; /**< \f$ \Omega_{0 b} \f$: baryons */
    double Omega0_ur; /**< \f$ \Omega_{0 \nu r} \f$: ultra-relativistic neutrinos */
    double Omega0_cdm;      /**< \f$ \Omega_{0 cdm} \f$: cold dark matter */
    double Omega0_chi;   /**< decaying darkmatter */
    double Omega0_lambda;    /**< \f$ \Omega_{0_\Lambda} \f$: cosmological constant */
    double Rhocrit0;
    /************************************/
    double Gamma0;      /* Gamma_0 for decaying dark matter in Mpc^{-1}*/
    double Gammad; /*Power index*/
    double DNeff; /*Delta N_eff at a= a_tr. Assumed to be constant(?)*/
    double atr;
};

double coeff[3] = {3.0/2.0, -2.0, 1.0/2.0};

/*Conformal Hubble constant in eV^-1*/
double H(double x, double rho_ddm, double rho_CFT, struct background *pba)
{
    double rhotot;
    rhotot = rho_CFT + rho_ddm + pba->Rhocrit0*((pba->Omega0_b + pba->Omega0_cdm)*exp(-3*x) + (pba->Omega0_g + pba->Omega0_ur)*exp(-4*x) + pba->Omega0_lambda*exp(0.039*x));
    if (rhotot <= 0){
        printf("Negative energy at a = %e!\n", exp(x));
        exit(1);
    }
    double log_H2 = 2*x + log(rhotot) + log(8*_PI_*_G_/3.0);
    return exp(0.5*log_H2)/(1e3*_c_);
}

/*Jacobian for BDF*/
void jacobian(double x, double *rho, double a, struct background *pba)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            pba->jac[i][j] = 0;
        }
    }
    double Hub = H(x, exp(rho[0]), exp(rho[1]), pba);
    double h = -_H_STEP_;
    double Gamma0 = pba->Gamma0;
    double Gammad = pba->Gammad;

    pba->jac[0][0] = h/a;
    pba->jac[0][1] = 0.0;
    pba->jac[1][0] = exp((Gammad+1)*x + rho[0])*pow(h,2)*Gamma0/(pow(a,2)*exp(rho[1])*Hub + a*exp((Gammad+1)*x+rho[0])*h*Gamma0);
    pba->jac[1][1] = exp(rho[1])*h*Hub/(a*exp(rho[1])*Hub + exp((Gammad+1)*x + rho[0])*h*Gamma0);
}

/*Evolution equation of rho*/
void deriv(double x, double *rho, double *drho, struct background *pba)
{
    double h = H(x, exp(rho[0]), exp(rho[1]), pba);
    drho[0] = -3 - (pba->Gamma0/h)*exp((pba->Gammad+1)*x);
    drho[1] = -4 + (pba->Gamma0/h)*exp((pba->Gammad+1)*x)*exp(rho[0] - rho[1]);
}

/*BDF2 solver*/
void bdf2(double x_next ,double y_prev[], double y_prev2[],
          double y_next[], int step, struct background *pba)
{
    double y_guess[2];
    double h = -_H_STEP_;
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double F[2];
        double df[2];
        double delta[2] = {0.0,0.0};
        deriv(x_next, y_guess, df, pba);
        jacobian(x_next, y_guess, coeff[0], pba);

        F[0] = coeff[0]*y_guess[0]/h + coeff[1]*y_prev[0]/h + coeff[2]*y_prev2[0]/h - df[0];
        F[1] = coeff[0]*y_guess[1]/h + coeff[1]*y_prev[1]/h + coeff[2]*y_prev2[1]/h - df[1];

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] += -F[j]*pba->jac[i][j];
            }
            y_guess[i] += delta[i];
        }
        
        double max_delta = _TOL_*1e-20;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i]) > max_delta) max_delta = fabs(delta[i]); 
        }
        if (max_delta < _TOL_) break;
    }
    pba->rho_chi_table[step] = y_guess[0];
    pba->rho_cft_table[step] = y_guess[1];
    pba->Gamma_table[step] = pba->Gamma0*exp((pba->Gammad + 1)*x_next);
    pba->H_table[step] = H(x_next, y_guess[0], y_guess[1], pba);
    memcpy(y_next, y_guess, sizeof(double) * 2);
}
/************************************************************************************* */

/*Functions for Interpolation(Hermite)*/
double dfdh(int idx, double *rho)
{
    return (rho[idx-1] - rho[idx+1])/(2.0*_H_STEP_);
}

double d2fdh2(int idx, double *rho)
{
    return (rho[idx-2] - 2*rho[idx-1] + rho[idx])/pow(_H_STEP_, 2);
}

double divdiff(int *idx, int n, double *rho, double *loga)
{
    if (n == 2){
        if (idx[1] == idx[0]) return dfdh(idx[0], rho);
        else return (rho[idx[1]] - rho[idx[0]])/(loga[idx[1]] - loga[idx[0]]);
    }
    if (n == 3 && (idx[0] == idx[1] && idx[1] == idx[2])) return d2fdh2(idx[0], rho)/2.0;
    else{
        int *subidx1 = (int *)malloc(sizeof(int) * (n-1));
        int *subidx2 = (int *)malloc(sizeof(int) * (n-1));
        double res;
        for (int i = 1; i < n; i++){
            subidx1[i-1] = idx[i];
        }
        for (int i = 0; i < n-1; i++){
            subidx2[i] = idx[i];
        }
        res = (divdiff(subidx1, n-1, rho, loga) - divdiff(subidx2, n-1, rho, loga))/(loga[idx[n-1]] - loga[idx[0]]);
        free(subidx1); free(subidx2);
        return res;
    }
}

/*input a is indeed log(a)*/
double Hermite(double a, double *rho, double *loga)
{
    int N = sizeof(loga)/sizeof(double);
    double res = 0.0;
    double **arr;
    arr = (double **)malloc(sizeof(double *) * 9);
    for (int i = 0; i < 9; i++){
        arr[i] = (double *)malloc(sizeof(double) * (9-i));
    }

    int idx1, idx2, idx3;
    double x1, x2, x3;
    double diff = 1e10;
    for (int i = 0; i < N; i++){
        if ((fabs(loga[i] - a)) < diff){
        idx2 = i;
        diff = fabs(loga[i] - a);
        }
    }
    if (idx2 == N-1) return rho[N-1];
    if (idx2 == 0) return rho[0];
    idx1 = idx2+1;
    idx3 = idx2-1;
    int group[9] = {idx1, idx1, idx1, idx2, idx2, idx2, idx3, idx3, idx3};
    x1 = loga[idx1]; x2 = loga[idx2]; x3 = loga[idx3];
    for (int i = 0; i < 9; i++){
        arr[0][i] = rho[group[i]];
    }

    for (int i = 1; i < 9; i++){
        for (int j = 0; j < 9-i; j++){
            int *idxs = (int *)malloc(sizeof(int) * (i+1));
            for (int k = 0; k < i+1; k++){
                idxs[k] = group[j+k];
            }
            arr[i][j] = divdiff(idxs, i+1, rho, loga);
            free(idxs);
        }
    }
    res += arr[0][0];
    for (int k = 1; k < 9; k++){
        double tmp = arr[k][0];
        for (int j = 0; j <= k-1; j++){
            tmp *= (a-loga[group[j]]);
        }
        res += tmp;
    }

    for (int i = 0; i < 9; i++){
        free(arr[i]);
    }
    free(arr);

    return res;
}

double interpolation(double *loga, double *rho, double a)
{
    if (a > exp(loga[0])){
        return rho[0] + 3.0 * log(exp(loga[0])/exp(a));
    }
    else if (a <= log(_AINIT_)){
        printf("a must be larger than ainit.\n");
        return 1e100;
    }
    else{
        return Hermite(a, rho, loga);
    }
}

int background_solve_my_component(struct background *pba) {
    pba->a_size = 100;
    class_alloc(pba->a_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_chi_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_cft_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->Gamma_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->H_table, pba->a_size * sizeof(double), pba->error_message);
  
    double a_ini = pba->atr;
    double a_end = _AINIT_;
    double h = -_H_STEP_;
    double RHOC = 3.0 * pow(pba->H0, 2.)/(8.0 * _PI_ * _G_ * pow(_c_, -2.));
    pba->Rhocrit0 = RHOC;
    double x = log(pba->atr);

    //fill a_table with log-spaced values
    for (int i = 0; i < pba->a_size; i++) pba->a_table[i] = log(a_ini) + i*h;

    double rho_cft_atr = RHOC*(pba->Omega0_g)*pow(pba->atr, -4.)*(7./8.)*pow(4./11., 4./3.)*(pba->DNeff);
    double rho_chi_atr = RHOC*(pba->Omega0_chi)*pow(pba->atr, -3.);
    double rho_prev2[2] = {log(rho_chi_atr) - 3*h, log(rho_cft_atr) - 4*h};
    double rho_prev[2];
    double rho_next[2];
    
    memcpy(rho_prev, rho_prev2, sizeof(double)*2);
    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double df[2];
        double F[2];
        double delta[2] = {0.0, 0.0};
        deriv(x+h, rho_prev, df, pba);
        jacobian(x+h, rho_prev, 1, pba);
        for (int i = 0; i < 2; i++){
            F[i] = rho_prev[i]/h - rho_prev2[i]/h - df[i];
        }
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] += -F[j]*pba->jac[i][j];
            }
            rho_prev[i] += delta[i];
        }

        double max_delta = _TOL_*1e-20;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i]) > max_delta) max_delta = fabs(delta[i]); 
        }
        if (max_delta < _TOL_) break;
    }

    /*Main BDF2 loop*/
    for (int i = 0; i < pba->a_size; i++){
        x += h;
        bdf2(x, rho_prev, rho_prev2, rho_next, i, pba);
        memcpy(rho_prev2, rho_prev, sizeof(double) * 2);
        memcpy(rho_prev, rho_next, sizeof(double) * 2);
    }
    return 0;
}

int main()
{
    double Gamma0 = 1e3; /*Gyr^-1*/
    double H0 = 70.0; /*km/s/Mpc*/

    struct background pba;
    /*Set background parameters*/
    pba.atr = 1e-3;
    pba.DNeff = 0.6;
    pba.Gamma0 = Gamma0 * _Gyr_over_Mpc_; /*(Gyr^-1) * Conversion factor = Mpc^-1*/
    pba.Gammad = -3.5;
    pba.H0 = H0/(1e3*_c_); /*(km/s/Mpc) * Conversion factor = Mpc^-1*/
    pba.Omega0_b = 0.02238280; /*Baryon*/
    pba.Omega0_cdm = 0.0; /*CDM*/
    pba.Omega0_chi = 0.1201075; /*Decaying dark matter*/
    pba.Omega0_g = 8e-5; /*Photon*/
    pba.Omega0_lambda = 0.67; /*Dark energy*/
    pba.Omega0_ur = 1e-6; /*Neutrinos*/

    double loga[1000];
    for (int i = 0; i < 1000; i++){
        loga[i] = log(pba.atr) + i*_H_STEP_;
    }

    if (background_solve_my_component(&pba) != _SUCCESS_) {
        printf("Error: %s\n", pba.error_message);
        return 1;
      }
    
    printf("Interpolated rho_x(a) at selected points:\n");
    FILE *output, *params;
    output = fopen("output.dat", "w");
    params = fopen("params.dat", "w");
    for (int i = 0; i < 1000; i++){
        fprintf(output, "%e\t%e\t%e\t%e\t%e\n",loga[i], interpolation(pba.a_table, pba.rho_chi_table, loga[i]), interpolation(pba.a_table, pba.rho_cft_table, loga[i]),
        interpolation(pba.a_table, pba.Gamma_table, loga[i]),interpolation(pba.a_table, pba.H_table, loga[i]));
    }

    
      // Free memory
    free(pba.a_table);
    free(pba.rho_chi_table);
    free(pba.rho_cft_table);
    free(pba.H_table);
    free(pba.Gamma_table);
    
    fprintf(params, "%s\t%e\n", "aeq", _AEQ_);
    fprintf(params, "%s\t%e\n", "Gammad", pba.Gammad);
    fprintf(params, "%s\t%e\n", "Gamma0", pba.Gamma0);
    fprintf(params, "%s\t%e\n", "DNeff", pba.DNeff);
    fprintf(params, "%s\t%e\n", "atr", pba.atr);
    fprintf(params, "%s\t%e\n", "ainit", _AINIT_);
    fclose(params); fclose(output);

    return 0;
}