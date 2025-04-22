#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

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
//#define _G_ 6.67428e-11    
#define _G_ 2.75e-115         /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */

#define _AEQ_ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define _PI_ 3.141592
#define _AINIT_ 1e-10

#define _RTOL_ 1e-6
#define _ATOL_ 1e-8
#define _MAX_ITER_ 20

double coeff[3] = {1.5, -2.0, 0.5};
double coeff2[4] = {11./6., -3., 1.5, -1./3.};
double coeff3[5] = {25./12., -4., 3., -4./3., 1./4.};

struct background {
    double *a_table;
    double *rho_chi_table;
    double *rho_cft_table;
    double *Gamma_table;
    double *H_table;
    double jac[2][2];
    int a_size;
    char error_message[1024];
    double steph;

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

/* Prototype for RHS function: 
 *   t   : current time
 *   y   : state vector of length neq
 *   dydt: output derivative vector of length neq
 */
typedef void (*RHSFunc)(double x, const double *y, double *dydt, struct background *pba ,int neq);

/*Hubble parameter in Mpc^-1*/
double H(double x, const double *rho, struct background *pba)
{
    double rhotot;
    rhotot = (exp(rho[0]) + exp(rho[1])) + pba->H0*pba->H0*((pba->Omega0_b + pba->Omega0_cdm)*exp(-3.*x) + (pba->Omega0_g + pba->Omega0_ur)*exp(-4.*x) + pba->Omega0_lambda);
    if (rhotot <= 0){
        printf("Negative energy at a = %e!\n", exp(x));
    }
    return sqrt(rhotot);
}

/*d(log rho)/d(log a)*/
void deriv(double x, const double *rho, double *drho, struct background *pba, int neq)
{
    double h = H(x, rho, pba);
    drho[0] = -3 - (pba->Gamma0/h)*exp(pba->Gammad*x);
    drho[1] = -4 + exp(rho[0] - rho[1])*(pba->Gamma0/h)*exp(pba->Gammad*x);
}

double jacdet(struct background *pba)
{
    return pba->jac[1][1]*pba->jac[0][0] - pba->jac[1][0]*pba->jac[0][1];
}

void jacobian(int n, double x, double *rho, struct background *pba)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            pba->jac[i][j] = 0.;
        }
    }
    double h;
    double fplus[2];
    double fminus[2];
    double dfplus[2];
    double dfminus[2];
    double rhoplus[2];
    double rhominus[2];
    double F[2];
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            memcpy(rhoplus, rho, sizeof(double)*2);
            memcpy(rhominus, rho, sizeof(double)*2);
            h = fmax(sqrt(DBL_EPSILON),_RTOL_)*(_ATOL_+fabs(rho[j]));
            rhoplus[j] += h;
            rhominus[j] -= h;
            deriv(x, rhoplus, dfplus, pba, 2);
            deriv(x, rhominus, dfminus, pba, 2);
            if (n == 2){
                fplus[i] = coeff[0]*rhoplus[i] - pba->steph*dfplus[i];
                fminus[i] = coeff[0]*rhominus[i] - pba->steph*dfminus[i];
            }
            else if (n == 3){
                fplus[i] = coeff2[0]*rhoplus[i] - pba->steph*dfplus[i];
                fminus[i] = coeff2[0]*rhominus[i] - pba->steph*dfminus[i];
            }
            else if (n == 1){
                fplus[i] = rhoplus[i] - pba->steph*dfplus[i];
                fminus[i] = rhominus[i] - pba->steph*dfminus[i];
            }
            else if (n == 4){
                fplus[i] = coeff3[0]*rhoplus[i] - pba->steph*dfplus[i];
                fminus[i] = coeff3[0]*rhominus[i] - pba->steph*dfminus[i];
            }
            pba->jac[i][j] = (fplus[i] - fminus[i])/(2.*h);
        }
    }
    double J = jacdet(pba);
    pba->jac[0][0] = pba->jac[1][1]/J;
    pba->jac[0][1] = -pba->jac[0][1]/J;
    pba->jac[1][0] = -pba->jac[1][0]/J;
    pba->jac[1][1] = pba->jac[0][0]/J;
}

static void get_ndf_coefs(int k, double *alpha, double *beta) {
    double theta;
    switch(k) {
        case 1:
            /* 1-step NDF == backward Euler */
            alpha[0] = 1.0; alpha[1] = -1.0;
            *beta = 1.0;
            break;
        case 2:
            theta = 2.0/3.0;
            /* alpha[0..2] */
            alpha[0] = 1.5 + theta;      /* 3/2 + 2/3 = 13/6 */
            alpha[1] = -2.0*(1.0+theta);  /* -2*(5/3) = -10/3 */
            alpha[2] = 0.5 + theta;       /* 1/2 + 2/3 = 7/6 */
            *beta    = 1.0;
            break;
        case 3:
            theta = 3.0/4.0;
            /* coefficients from Shampine & Reichelt */
            alpha[0] = 11.0/6.0 + theta;  /* approx 11/6 + 3/4 = 2.5833 */
            alpha[1] = -3.0*(1.0+theta);  /* -3*(7/4) = -5.25 */
            alpha[2] = 1.5 + 3*theta/2.0;  /* 3/2 + 9/8 = 2.625 */
            alpha[3] = -(1.0/3.0 + theta);/* - (1/3+3/4)= -1.0833 */
            *beta    = 1.0;
            break;
        default:
            fprintf(stderr, "Error: Unsupported NDF order %d\n", k);
            exit(EXIT_FAILURE);
    }
}

/* Solve ODE dy/dt = f(t,y) from t0 to tf with fixed step h and order k (1<=k<=3)
 * y0: initial state vector length neq
 * sol: output array of size (nsteps+1)*neq, row-major [step][i]
 * Returns number of steps taken
 */
void ndf_solve(RHSFunc f, int neq,
              double t0, double tf, double h,
              int k, const double *y0, double *tout, struct background *pba) {
    int max_history = k+1;
    double **Yhist = malloc(max_history * sizeof(double*));
    for(int i = 0; i < max_history; i++)
        Yhist[i] = malloc(neq * sizeof(double));
    double *dydt = malloc(neq * sizeof(double));
    double *alpha = malloc((k+1) * sizeof(double));
    double beta;

    /* initialize history with initial condition */
    for(int i = 0; i < neq; i++) Yhist[0][i] = y0[i];
    double t = t0;

    /* simple forward initialization for Yhist[1..k-1] */
    for(int j = 1; j < k; j++) {
        f(t, Yhist[j-1], dydt, pba, neq);
        /* step with same h sign for consistency */
        t += h;
        for(int i = 0; i < neq; i++)
            Yhist[j][i] = Yhist[j-1][i] + h * dydt[i];
    }

    int max_steps = pba->a_size;
    for(int step = 0; step < max_steps; step++) {
        /* adjust last step size to land exactly on tf */
        if((h > 0 && t + h > tf) || (h < 0 && t + h < tf)) {
            h = tf - t;
        }
        get_ndf_coefs(k, alpha, &beta);
        double *ynew = malloc(neq * sizeof(double));
        memcpy(ynew, Yhist[0], neq*sizeof(double));

        for (int iter = 0; iter < _MAX_ITER_; iter++){
            f(t, ynew, dydt, pba, neq);
            double *F = malloc(neq * sizeof(double));
            double *delta = malloc(neq * sizeof(double));
            delta[0] = 0.0;
            delta[1] = 0.0;
            jacobian(k, t, ynew, pba);
            F[0] = -h*beta*dydt[0];
            F[1] = -h*beta*dydt[1];
            for(int i = 0; i < neq; i++) {
                for(int j = 0; j <= k; j++) F[i] += alpha[j] * Yhist[j][i];
            }
            for (int i = 0; i < neq; i++){
                for (int j = 0; j < neq; j++){
                    delta[i] -= pba->jac[i][j]*F[j];
                }
                ynew[i] += delta[i];
            }
            double max_delta = 0.0;
            for (int i = 0; i < 2; i++){
                if (fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(ynew[i] - delta[i])) > max_delta) max_delta = fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(ynew[i] - delta[i])); 
            }
            free(F); free(delta);
            if (max_delta < 1) break;
        }

        /* shift history */
        for(int j = k; j > 0; j--) memcpy(Yhist[j], Yhist[j-1], neq * sizeof(double));
        memcpy(Yhist[0], ynew, neq * sizeof(double));
        free(ynew);

        t += h;
        /* store output */
        pba->rho_chi_table[step] = Yhist[0][0];
        pba->rho_cft_table[step] = Yhist[0][1];
        pba->a_table[step] = t;
        pba->Gamma_table[step] = pba->Gamma0*pow(exp(t),pba->Gammad);
        pba->H_table[step] = H(t, Yhist[0], pba);

        /* break if reached or passed tf */
        if((h > 0 && t >= tf) || (h < 0 && t <= tf)) break;
    }

    /* cleanup */
    for(int i = 0; i < max_history; i++) free(Yhist[i]);
    free(Yhist); free(dydt); free(alpha);
}

/*Functions for Interpolation(Hermite)*/
double dfdh(int idx, double *rho, struct background *pba)
{
    return (rho[idx-1] - rho[idx+1])/(2.0*pba->steph);
}

double d2fdh2(int idx, double *rho, struct background *pba)
{
    return (rho[idx-2] - 2*rho[idx-1] + rho[idx])/(pba->steph * pba->steph);
}

double divdiff(int *idx, int n, double *rho, struct background *pba)
{
    if (n == 2){
        if (idx[1] == idx[0]) return dfdh(idx[0], rho, pba);
        else return (rho[idx[1]] - rho[idx[0]])/(pba->a_table[idx[1]] - pba->a_table[idx[0]]);
    }
    if (n == 3 && (idx[0] == idx[1] && idx[1] == idx[2])) return d2fdh2(idx[0], rho, pba)/2.0;
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
        res = (divdiff(subidx1, n-1, rho, pba) - divdiff(subidx2, n-1, rho, pba))/(pba->a_table[idx[n-1]] - pba->a_table[idx[0]]);
        free(subidx1); free(subidx2);
        return res;
    }
}

/*input a is indeed log(a)*/
double Hermite(struct background *pba, double *rho, double a)
{
    int N = pba->a_size;
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
        if ((fabs(pba->a_table[i] - a)) < diff){
        idx2 = i;
        diff = fabs(pba->a_table[i] - a);
        }
    }
    if (idx2 == N-1) return rho[N-1];
    if (idx2 == 0) return rho[0];
    idx1 = idx2+1;
    idx3 = idx2-1;
    int group[9] = {idx1, idx1, idx1, idx2, idx2, idx2, idx3, idx3, idx3};
    x1 = pba->a_table[idx1]; x2 = pba->a_table[idx2]; x3 = pba->a_table[idx3];
    for (int i = 0; i < 9; i++){
        arr[0][i] = rho[group[i]];
    }

    for (int i = 1; i < 9; i++){
        for (int j = 0; j < 9-i; j++){
            int *idxs = (int *)malloc(sizeof(int) * (i+1));
            for (int k = 0; k < i+1; k++){
                idxs[k] = group[j+k];
            }
            arr[i][j] = divdiff(idxs, i+1, rho, pba);
            free(idxs);
        }
    }
    res += arr[0][0];
    for (int k = 1; k < 9; k++){
        double tmp = arr[k][0];
        for (int j = 0; j <= k-1; j++){
            tmp *= (a-pba->a_table[group[j]]);
        }
        res += tmp;
    }

    for (int i = 0; i < 9; i++){
        free(arr[i]);
    }
    free(arr);

    return res;
}
/*Interpolation Done*/

double interpolation(struct background *pba, double *rho, double a)
{
    if (a > log(pba->atr)){
        return log(exp(rho[0])*pow(pba->atr/exp(a), 3));
    }
    else if (a < log(_AINIT_)){
        printf("a must be larger than ainit.\n");
        return 1e100;
    }
    else{
        return Hermite(pba, rho, a);
    }
}

double GammaChi(double x, struct background *pba)
{
    if (x <= log(pba->atr)) return pba->Gamma0*pow(exp(x)/pba->atr, pba->Gammad);
    else return 0.0;
}

int background_solve_my_component(struct background *pba) {
    pba->a_size = 10000;
    class_alloc(pba->a_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_chi_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_cft_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->Gamma_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->H_table, pba->a_size * sizeof(double), pba->error_message);
  
    double a_ini = pba->atr;
    double a_end = _AINIT_;
    pba->steph = -((log(a_ini) - log(a_end))/(double)(pba->a_size-1));
    double x = log(a_ini);

    //Here we are using normalization where 8*pi*G/3 = 1. 
    double rho_cft_atr = pba->H0*pba->H0*(pba->Omega0_g)*pow(pba->atr, -4.)*(7./8.)*pow(4./11., 4./3.)*(pba->DNeff);
    double rho_chi_atr = pba->H0*pba->H0*(pba->Omega0_chi)*pow(pba->atr, -3.);
    double rho0[2] = {log(rho_chi_atr), log(rho_cft_atr)};

    ndf_solve(deriv, 2, x, log(a_end), pba->steph, 3, rho0, pba->a_table, pba);
    return 0;
}

int main()
{
    double Gamma0_input = 1e3; /*km/s/Mpc*/
    double H0 = 70.0; /*km/s/Mpc. This is test input. Real code will get theta as input.*/

    struct background pba;
    /*Set background parameters*/
    pba.atr = 1e-3;
    pba.DNeff = 0.5;
    pba.Gamma0 = Gamma0_input*1e3/_c_; /*(km/s/Mpc) * Conversion factor = Mpc^-1*/
    pba.Gammad = 0.;
    pba.H0 = H0*1e3/_c_; /*(km/s/Mpc) * Conversion factor = Mpc^-1*/
    pba.Omega0_b = 0.02238280; /*Baryon*/
    pba.Omega0_cdm = 0.1; /*CDM*/
    pba.Omega0_chi = 0.3; /*Decaying dark matter*/
    pba.Omega0_g = 5e-5; /*Photon*/
    pba.Omega0_lambda = 0.67; /*Dark energy*/
    pba.Omega0_ur = 1e-6; /*Neutrinos*/

    if (background_solve_my_component(&pba) != _SUCCESS_) {
        printf("Error: %s\n", pba.error_message);
        return 1;
      }
    
    printf("Interpolated rho_x(a) at selected points:\n");
    FILE *output, *params;
    output = fopen("output.dat", "w");
    params = fopen("params.dat", "w");

    int N = 10000;
    double *loga = (double *)malloc(sizeof(double)*N);
    double ai = _AINIT_;
    double aend = 1.0;
    double h = (log(aend) - log(ai))/N;
    for (int i = 0; i < N; i++){
        loga[i] = log(_AINIT_) + i*h;
    }
    
    
    //Interpolation from computed samples
    /*
    for (int i = 0; i < N; i++){
        double rho[2];
        double rho_chi = interpolation(&pba, pba.rho_chi_table, loga[i]);
        double rho_cft = interpolation(&pba, pba.rho_cft_table, loga[i]);
        rho[0] = rho_chi; rho[1] = rho_cft;
        fprintf(output, "%e\t%e\t%e\t%e\t%e\n",loga[i], rho_chi, rho_cft, GammaChi(loga[i], &pba),H(loga[i], rho, &pba));
    }*/
    printf("Done.\n");
    
    for (int i = 0; i < pba.a_size; i++){
        fprintf(output, "%.20g\t%.20g\t%.20g\t%.20g\t%.20g\n", pba.a_table[i], pba.rho_chi_table[i], pba.rho_cft_table[i], pba.Gamma_table[i], pba.H_table[i]);
    }

    
      // Free memory
    free(pba.a_table);
    free(pba.rho_chi_table);
    free(pba.rho_cft_table);
    free(pba.H_table);
    free(pba.Gamma_table);
    free(loga);
    
    fprintf(params, "%s\t%e\n", "aeq", _AEQ_);
    fprintf(params, "%s\t%e\n", "Gammad", pba.Gammad);
    fprintf(params, "%s\t%e\n", "Gamma0", Gamma0_input);
    fprintf(params, "%s\t%e\n", "DNeff", pba.DNeff);
    fprintf(params, "%s\t%e\n", "atr", pba.atr);
    fprintf(params, "%s\t%e\n", "ainit", _AINIT_);
    fclose(params); fclose(output);

    return 0;
}