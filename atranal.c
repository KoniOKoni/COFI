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
                        
#define _c_ 2.99792458e8            /**< c in m/s */
//#define _G_ 6.67428e-11    
#define _G_ 2.75e-115         /**< Newton constant in m^3/Kg/s^2 */
#define _eV_ 1.602176487e-19        /**< 1 eV expressed in J */
#define _year_to_Mpc_ 3.07e-7

#define _AEQ_ 2.93e-4 /*Scale factor at matter-radiation equality*/
#define _PI_ 3.141592
#define _AINIT_ 1e-10

#define _RTOL_ 1e-6
#define _ATOL_ 1e-8
#define _MAX_ITER_ 100

#define _TRUE_ 1
#define _FALSE_ 0

//Backgorund structure
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
    int tablesize;

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
    double DDM_decay_time; //Decay time 'tau' of DDM in years.
    double f_DDM;

    double rho_chi_init;
    double rho_cft_init;
};

#define BAR_WIDTH 50

void print_progress(double fraction) {
    int filled = (int)(fraction * BAR_WIDTH);
    int empty  = BAR_WIDTH - filled;
    int percent = (int)(fraction * 100);

    printf("\r%3d%% [", percent);
    for (int i = 0; i < filled; i++) putchar('#');
    for (int i = 0; i < empty; i++) putchar(' ');
    putchar(']');

    fflush(stdout);
}

//Integrator coefficients
double coeff[3] = {1.5, -2.0, 0.5};
double coeff2[4] = {11./6., -3., 1.5, -1./3.};
double coeff3[5] = {25./12., -4., 3., -4./3., 1./4.};

//Determinant of the Jacobian
double jacdet(struct background *pba)
{
    return pba->jac[1][1]*pba->jac[0][0] - pba->jac[1][0]*pba->jac[0][1];
}

/*Hubble parameter in Mpc^-1*/
double H(double x, double *rho, struct background *pba)
{
    double rhotot;
    rhotot = rho[0] + rho[1] + pow(pba->H0, 2)*(pba->Omega0_b*exp(-3.*x) + (pba->Omega0_g + pba->Omega0_ur)*exp(-4.*x) + pba->Omega0_lambda);
    return sqrt(rhotot);
}

/*d rho/d(log a)*/
void deriv(double x, double *rho, double *drho, struct background *pba)
{
    double h = H(x, rho, pba);
    if (x < log(pba->atr)){
        drho[0] = -3*rho[0] - (pba->Gamma0/h)*rho[0]*exp(x*pba->Gammad);
        drho[1] = -4*rho[1] + (pba->Gamma0/h)*rho[0]*exp(x*pba->Gammad);
    }
    else{
        drho[0] = -3*rho[0];
        drho[1] = -4*rho[1];
    }
}

//Compute Jacobian numerically.
//This function actually gives J^{-1} which is finally needed.
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
            h = sqrt(DBL_EPSILON)*(1+fabs(rho[j]));
            rhoplus[j] += h;
            rhominus[j] -= h;
            deriv(x, rhoplus, dfplus, pba);
            deriv(x, rhominus, dfminus, pba);
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
    double tmp[2][2];

    tmp[0][0] = pba->jac[1][1]/J;
    tmp[1][1] = pba->jac[0][0]/J;
    tmp[0][1] = -pba->jac[1][0]/J;
    tmp[1][0] = -pba->jac[0][1]/J;

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            pba->jac[i][j] = tmp[i][j];
        }
    }
}

//BDF solvers
/* ****************************************************************************************** */
void bdf2(double x_next ,double y_prev[], double y_prev2[],
          double y_next[], int step, struct background *pba)
{
    double y_guess[2];
    double h = pba->steph;
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double F[2];
        double df[2];
        double delta[2] = {0.0,0.0};
        deriv(x_next, y_guess, df, pba);
        jacobian(2, x_next, y_guess, pba);

        F[0] = coeff[0]*y_guess[0] + coeff[1]*y_prev[0] + coeff[2]*y_prev2[0] - h*df[0];
        F[1] = coeff[0]*y_guess[1] + coeff[1]*y_prev[1] + coeff[2]*y_prev2[1] - h*df[1];
        /*for (int i = 0; i < 2; i++){
            printf("%e\t%e\t%e\n",y_guess[i], y_prev[i], y_prev2[i]);
        }*/

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] -= F[j]*pba->jac[i][j];
            }
            y_guess[i] += delta[i];
        }
        
        double max_delta = 0.0;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])) > max_delta) max_delta = fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])); 
        }
        if (max_delta < 1) break;
    }
    pba->rho_chi_table[step] = fmax(y_guess[0],0.);
    pba->rho_cft_table[step] = fmax(y_guess[1],0.);
    pba->Gamma_table[step] = pba->Gamma0*pow(exp(x_next), pba->Gammad);
    pba->H_table[step] = H(x_next, y_guess, pba);
    pba->tablesize = step;
    memcpy(y_next, y_guess, sizeof(double) * 2);
}

void bdf3(double x_next ,double y_prev[], double y_prev2[], double y_prev3[],
          double y_next[], int step, struct background *pba)
{
    double y_guess[2];
    double h = pba->steph;
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double F[2];
        double df[2];
        double delta[2] = {0.0,0.0};
        deriv(x_next, y_guess, df, pba);
        jacobian(3, x_next, y_guess, pba);
        //printf("a = %e\tdrho(chi) = %e\tdrho(cft) = %e\n",exp(x_next), df[0], df[1]); 

        for (int i = 0; i < 2; i++){
            F[i] = y_guess[i]*coeff2[0] + y_prev[i]*coeff2[1] + y_prev2[i]*coeff2[2] + y_prev3[i]*coeff2[3] - h*df[i];
        }

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] -= F[j]*pba->jac[i][j];
            }
            y_guess[i] += delta[i];
        }
        
        double max_delta = 0.0;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])) > max_delta) max_delta = fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])); 
        }
        if (max_delta < 1) break;
    }
    pba->rho_chi_table[step] = fmax(y_guess[0], 0.);
    pba->rho_cft_table[step] = fmax(y_guess[1], 0.);
    pba->Gamma_table[step] = pba->Gamma0*pow(exp(x_next), pba->Gammad);
    pba->H_table[step] = H(x_next, y_guess, pba);
    pba->tablesize = step;
    memcpy(y_next, y_guess, sizeof(double) * 2);
}

void bdf4(double x_next ,double y_prev[], double y_prev2[], double y_prev3[], double y_prev4[],
          double y_next[], int step, struct background *pba)
{
    double y_guess[2];
    double h = pba->steph;
    memcpy(y_guess, y_prev, sizeof(double)*2);

    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double F[2];
        double df[2];
        double delta[2] = {0.0,0.0};
        deriv(x_next, y_guess, df, pba);
        jacobian(4, x_next, y_guess, pba);
        //printf("a = %e\tdrho(chi) = %e\tdrho(cft) = %e\n",exp(x_next), df[0], df[1]); 

        for (int i = 0; i < 2; i++){
            F[i] = y_guess[i]*coeff3[0] + y_prev[i]*coeff3[1] + y_prev2[i]*coeff3[2] + y_prev3[i]*coeff3[3] + y_prev4[i]*coeff3[4] - h*df[i];
        }

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] -= F[j]*pba->jac[i][j];
            }
            y_guess[i] += delta[i];
        }
        
        double max_delta = 0.0;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])) > max_delta) max_delta = fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(y_guess[i])); 
        }
        if (max_delta < 1) break;
    }
    //if (y_guess[0] < 0) y_guess[0] = 1e-30;
    //if (y_guess[1] < 0) y_guess[1] = 1e-30;
    pba->rho_chi_table[step] = fmax(y_guess[0], 0.);
    pba->rho_cft_table[step] = fmax(y_guess[1], 0.);
    pba->Gamma_table[step] = pba->Gamma0*pow(exp(x_next), pba->Gammad);
    pba->H_table[step] = H(x_next, y_guess, pba);
    pba->tablesize = step;
    //printf("a = %e\tdrho(chi) = %e\tdrho(cft) = %e\n",exp(x_next), df[0], df[1]); 
    memcpy(y_next, y_guess, sizeof(double) * 2);
}
/* ****************************************************************************************** */

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
    int N = pba->tablesize;
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

double GammaChi(double x, struct background *pba)
{
    if (x <= log(pba->atr)) return pba->Gamma0*pow(exp(x),pba->Gammad);
    else return 0.0;
}

int background_solve_my_component(struct background *pba) {
    pba->a_size = 1000;
    class_alloc(pba->a_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_chi_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->rho_cft_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->Gamma_table, pba->a_size * sizeof(double), pba->error_message);
    class_alloc(pba->H_table, pba->a_size * sizeof(double), pba->error_message);
  
    double a_ini = _AINIT_;
    double a_end = 1.0;
    pba->steph = ((log(a_end) - log(a_ini))/(double)(pba->a_size));
    double x = log(a_ini);

    //fill a_table with log-spaced values
    for (int i = 0; i < pba->a_size; i++) pba->a_table[i] = x + i*(pba->steph);

    //Here we are using normalization where 8*pi*G/3 = 1. Energy density is in unit of Mpc^-2.
    pba->rho_chi_init = pow(pba->H0, 2)*(pba->Omega0_cdm)*exp(-3*x)*(1. + pba->f_DDM);
    pba->rho_cft_init = 0.;

    //Arrays for BDF iteration.
    double rho_prev4[2] = {pba->rho_chi_init, pba->rho_cft_init};//Initial condition
    double rho_prev3[2];
    double rho_prev2[2];
    double rho_prev[2];
    double rho_next[2];

    //Save initial values
    pba->rho_chi_table[0] = rho_prev4[0]; pba->rho_cft_table[0] = rho_prev4[1]; pba->Gamma_table[0] = pba->Gamma0*pow(pba->atr, pba->Gammad); pba->H_table[0] = H(x, rho_prev4, pba);
    memcpy(rho_prev3, rho_prev4, sizeof(double)*2);

    //Newton's method once
    for (int iter = 0; iter < _MAX_ITER_; iter++){
        double df[2];
        double F[2];
        double delta[2] = {0.0, 0.0};
        deriv(x+pba->steph, rho_prev3, df, pba);
        jacobian(1, x+pba->steph, rho_prev3, pba);
        for (int i = 0; i < 2; i++){
            F[i] = rho_prev3[i] - rho_prev4[i] - pba->steph*df[i];
        }
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                delta[i] -= F[j]*pba->jac[i][j];
            }
            rho_prev3[i] += delta[i];
        }

        double max_delta = 0.0;
        for (int i = 0; i < 2; i++){
            if (fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(rho_prev3[i])) > max_delta) max_delta = fabs(delta[i])/(_ATOL_ + _RTOL_*fabs(rho_prev3[i])); 
        }
        if (max_delta < 1e-1) break;
    }
    x += pba->steph;
    pba->rho_chi_table[1] = rho_prev3[0]; pba->rho_cft_table[1] = rho_prev3[1]; pba->Gamma_table[1] = pba->Gamma0*pow(exp(x), pba->Gammad); pba->H_table[1] = H(x, rho_prev3, pba);
    x += pba->steph;
    bdf2(x, rho_prev3, rho_prev4, rho_prev2, 2, pba); //BDF2 once
    x += pba->steph;
    bdf3(x, rho_prev2, rho_prev3, rho_prev4, rho_prev, 3, pba); //BDF3 once
    /*Main BDF4 loop*/
    for (int i = 4; i < pba->a_size; i++){
        x += pba->steph;
        bdf4(x, rho_prev, rho_prev2, rho_prev3, rho_prev4, rho_next, i, pba);
        memcpy(rho_prev4, rho_prev3, sizeof(double) * 2);
        memcpy(rho_prev3, rho_prev2, sizeof(double) * 2);
        memcpy(rho_prev2, rho_prev, sizeof(double) * 2);
        memcpy(rho_prev, rho_next, sizeof(double) * 2);
    }
    return 1;
}

int main()
{
    double H0_input = 70.0; /*km/s/Mpc. This is test input. Real code will get theta as input.*/

    struct background pba;
    /*Set background parameters*/
    pba.Gammad = 0.;
    pba.Omega0_b = 0.02238280; /*Baryon*/
    pba.Omega0_cdm = 0.3; /*CDM*/
    pba.Omega0_g = 4.6e-5; /*Photon*/
    pba.Omega0_lambda = 0.67; /*Dark energy*/
    pba.Omega0_ur = 0.001; /*Neutrinos*/
    pba.H0 = H0_input*1e3/_c_; /*(km/s/Mpc) * Conversion factor = Mpc^-1*/
    pba.Gamma0 = 1./(pba.DDM_decay_time*_year_to_Mpc_); /*Mpc^-1*/

    FILE *rhoDMtau, *rhoDMatr;
    rhoDMtau = fopen("rhoDMtau.dat", "w");
    rhoDMatr = fopen("rhoDMatr.dat", "w");

    int N = 100, cnt = 0;
    int SUC;
    double DNeff, rho_cft_atr;

    double *f_table = malloc(sizeof(double) * N);
    double *tau_table = malloc(sizeof(double) * N);
    double *atr_table = malloc(sizeof(double) * N);

    for (int i = 0; i < N; i++){
        f_table[i] = log(1e-4) + i * (log(1e-1) - log(1e-4))/((double)N);
        tau_table[i] = log(1) + i * (log(100e9) - log(1))/((double)N);
        atr_table[i] = log(_AINIT_ * 10) + i * (log(1.) - log(_AINIT_*10))/((double)N);
    }

    pba.atr = 1./(1. + 1100.);
    //Fix a_tr and vary tau and rho_DM_init
    printf("Start rhoDM vs tau scanning....\n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            print_progress((double)cnt / (N*N));
            pba.f_DDM = exp(f_table[i]);
            pba.DDM_decay_time = exp(tau_table[j]);
            pba.Gamma0 = 1./(pba.DDM_decay_time*_year_to_Mpc_);
            SUC = background_solve_my_component(&pba);
            rho_cft_atr = Hermite(&pba, pba.rho_cft_table, log(pba.atr));
            if (rho_cft_atr > 0){
                DNeff = 11.*pow(22., 1/3)*rho_cft_atr*pow(pba.atr, 4)/(7.*pba.H0*pba.H0*pba.Omega0_g);
                if (DNeff >= 0.1 && DNeff <= 1.){
                    fprintf(rhoDMtau, "%.20g\t%.20g\t%.20g\n", pba.f_DDM, pba.DDM_decay_time,DNeff);
                }
            }
            cnt++;
            free(pba.a_table);
            free(pba.rho_chi_table);
            free(pba.rho_cft_table);
            free(pba.H_table);
            free(pba.Gamma_table);
        }
    }

    printf("Done.\n");
    printf("Start rhoDM vs atr scanning....\n");
    cnt = 0;
    pba.DDM_decay_time = 378000.;
    pba.Gamma0 = 1./(pba.DDM_decay_time*_year_to_Mpc_);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            print_progress((double)cnt / (N*N));
            pba.f_DDM = exp(f_table[i]);
            pba.atr = exp(atr_table[j]);
            SUC = background_solve_my_component(&pba);
            rho_cft_atr = Hermite(&pba, pba.rho_cft_table, log(pba.atr));
            if (rho_cft_atr > 0){
                DNeff = 11.*pow(22., 1/3)*rho_cft_atr*pow(pba.atr, 4)/(7.*pba.H0*pba.H0*pba.Omega0_g);
                if (DNeff >= 0.1 && DNeff <= 1.){
                    fprintf(rhoDMatr, "%.20g\t%.20g\t%.20g\n", pba.f_DDM, pba.atr,DNeff);
                }
            }
            cnt++;
            free(pba.a_table);
            free(pba.rho_chi_table);
            free(pba.rho_cft_table);
            free(pba.H_table);
            free(pba.Gamma_table);
        }
    }
    printf("Done.\n");
    fclose(rhoDMatr); fclose(rhoDMtau);
    free(f_table); free(tau_table); free(atr_table);

    return 0;
}