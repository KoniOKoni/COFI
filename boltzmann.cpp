#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <fstream>
#include <stdexcept>
#include <iomanip>

// Physical and conversion constants
static constexpr long double Mpc_over_m   = 3.085677581282e22;     // m per Mpc
static constexpr long double Gyr_over_Mpc = 3.06601394e2;          // Gyr per Mpc (c=1, Julian year)
static constexpr long double Conversion   = 980.392;               // conversion factor Gyr^-1 -> Mpc^-1
static constexpr long double c_light      = 2.99792458e8;          // m/s
static constexpr long double G_Newton    = 2.75e-115;              // m^3/kg/s^2
static constexpr long double eV_J        = 1.602176487e-19;        // J per eV
static constexpr long double AEQ         = 2.93e-4;                // scale factor at equality
static constexpr long double PI          = 3.14159265358979323846;
static constexpr long double AINIT       = 1e-8;
static constexpr long double RTOL        = 1e-6;
static constexpr long double ATOL        = 1e-8;
static constexpr int    MAX_ITER    = 100;

static const std::array<long double,3> coeff2_1 = {1.5, -2.0, 0.5};
static const std::array<long double,4> coeff3_1 = {11.0/6.0, -3.0, 1.5, -1.0/3.0};
static const std::array<long double,5> coeff4_1 = {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0};

class Background {
public:
    int a_size;
    long double steph;
    
    // Cosmological parameters
    long double H0, Omega0_g, Omega0_b, Omega0_ur, Omega0_cdm;
    long double Omega0_chi, Omega0_lambda, Rhocrit0;
    long double Gamma0, Gammad, DNeff, atr;

    // Data tables
    std::vector<long double> a_table, rho_chi_table, rho_cft_table, Gamma_table, H_table;
    std::array<std::array<long double,2>,2> jac;

    Background() = default;
    
    /** Solve background component; returns 0 on success */
    int solve() {
        a_size = 100000;
        a_table.resize(a_size);
        rho_chi_table.resize(a_size);
        rho_cft_table.resize(a_size);
        Gamma_table.resize(a_size);
        H_table.resize(a_size);

        // Step in log(a)
        steph = -(std::log(atr) - std::log(AINIT)) / static_cast<long double>(a_size);
        long double x = std::log(atr);
        for (int i = 0; i < a_size; ++i) {
            a_table[i] = x + i * steph;
        }

        // Initial densities at a_tr
        long double rho_cft_atr = H0 * H0 * Omega0_g * (7.0/8.0) * std::pow(4.0/11.0, 4.0/3.0) * DNeff;
        long double rho_chi_atr = H0 * H0 * Omega0_chi;

        // Reserve arrays for multi-step
        std::array<long double,2> rho_prev4 = {rho_chi_atr, rho_cft_atr};
        std::array<long double,2> rho_prev3;
        std::array<long double,2> rho_prev2;
        std::array<long double,2> rho_prev;
        std::array<long double,2> rho_next;

        // Save step 0
        rho_chi_table[0] = rho_prev4[0];
        rho_cft_table[0] = rho_prev4[1];
        Gamma_table[0] = Gamma0 * std::pow(atr, Gammad);
        H_table[0]     = H(x, rho_prev4);
        rho_prev3 = rho_prev4;

        // Single Newton step
        x += steph;
        newtonStep(x, rho_prev3, rho_prev4);
        rho_chi_table[1] = rho_prev3[0];
        rho_cft_table[1] = rho_prev3[1];
        Gamma_table[1]   = Gamma0 * std::pow(std::exp(x), Gammad);
        H_table[1]       = H(x, rho_prev3);

        // BDF2 first call
        x += steph;
        bdf2(x, rho_prev3, rho_prev4, rho_prev2, 2);
        rho_chi_table[2] = rho_prev2[0];
        rho_cft_table[2] = rho_prev2[1];
        Gamma_table[2]   = Gamma0 * std::pow(std::exp(x), Gammad);
        H_table[2]       = H(x, rho_prev2);

        // BDF3 first call
        x += steph;
        bdf3(x, rho_prev2, rho_prev3, rho_prev4, rho_prev, 3);
        rho_chi_table[3] = rho_prev[0];
        rho_cft_table[3] = rho_prev[1];
        Gamma_table[3]   = Gamma0 * std::pow(std::exp(x), Gammad);
        H_table[3]       = H(x, rho_prev);

        // Main BDF4 loop
        for (int i = 4; i < a_size; ++i) {
            x += steph;
            bdf4(x, rho_prev, rho_prev2, rho_prev3, rho_prev4, rho_next, i);
            rho_prev4 = rho_prev3;
            rho_prev3 = rho_prev2;
            rho_prev2 = rho_prev;
            rho_prev  = rho_next;
        }
        return 0;
    }

private:
    long double H(long double x, const std::array<long double,2>& rho) const {
        long double rhotot = rho[0]/std::exp(3.0*x)
                       + rho[1]/std::exp(4.0*x)
                       + H0*H0 * ((Omega0_b+Omega0_cdm)*std::exp(-3.0*x)
                                + (Omega0_g+Omega0_ur)*std::exp(-4.0*x)
                                + Omega0_lambda);
        if (rhotot <= 0) {
            std::cerr << "Negative energy at a = " << std::exp(x) << std::endl;
        }
        return std::sqrt(rhotot);
    }

    void deriv(long double x, const std::array<long double,2>& rho,
               std::array<long double,2>& drho) const {
        long double h_val = H(x, rho);
        drho[0] = -(Gamma0/h_val)*rho[0];
        drho[1] =  (Gamma0/h_val)*rho[0]*std::exp(x);
    }

    long double jacdet() const {
        return jac[1][1]*jac[0][0] - jac[1][0]*jac[0][1];
    }

    void computeJacobian(int n, long double x, const std::array<long double,2>& rho) {
        long double eps = std::numeric_limits<long double>::epsilon();
        std::array<long double,2> fplus, fminus, dfplus, dfminus;
        for (int i = 0; i < 2; ++i) jac[i].fill(0.0);

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                std::array<long double,2> rho_p = rho, rho_m = rho;
                long double h_step = eps;
                rho_p[j] += h_step;
                rho_m[j] -= h_step;
                deriv(x, rho_p, dfplus);
                deriv(x, rho_m, dfminus);

                auto fillF = [&](std::array<long double,2>& F, const std::array<long double,2>& r, const std::array<long double,2>& df)->void {
                    switch(n) {
                        case 1:
                            F[i] = r[i] - steph * df[i];
                            break;
                        case 2:
                            F[i] = coeff2_1[0]*r[i] - steph * df[i];
                            break;
                        case 3:
                            F[i] = coeff3_1[0]*r[i] - steph * df[i];
                            break;
                        case 4:
                            F[i] = coeff4_1[0]*r[i] - steph * df[i];
                            break;
                    }
                };
                fillF(fplus, rho_p, dfplus);
                fillF(fminus, rho_m, dfminus);
                jac[i][j] = (fplus[i] - fminus[i]) / (2.0 * h_step);
            }
        }
        // Invert 2x2 jacobian
        long double J = jacdet();
        std::array<std::array<long double,2>,2> inv;
        inv[0][0] =  jac[1][1]/J;
        inv[0][1] = -jac[0][1]/J;
        inv[1][0] = -jac[1][0]/J;
        inv[1][1] =  jac[0][0]/J;
        jac = inv;
    }

    void newtonStep(long double x, std::array<long double,2>& y, const std::array<long double,2>& y0) {
        for (int iter = 0; iter < MAX_ITER; ++iter) {
            std::array<long double,2> df, F, delta = {0.0, 0.0};
            deriv(x, y, df);
            computeJacobian(1, x, y);
            for (int i = 0; i < 2; ++i) {
                F[i] = y[i] - y0[i] - steph * df[i];
            }
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) delta[i] -= F[j] * jac[i][j];
                y[i] += delta[i];
            }
            long double maxd = 0;
            for (int i = 0; i < 2; ++i) {
                long double norm = std::fabs(delta[i]) / (ATOL + RTOL * std::fabs(y[i] - delta[i]));
                maxd = std::max(maxd, norm);
            }
            if (maxd < 1e-1) break;
        }
    }

    void bdf2(long double x_next,
              const std::array<long double,2>& y_prev,
              const std::array<long double,2>& y_prev2,
              std::array<long double,2>& y_next,
              int step) {
        std::array<long double,2> y_guess = y_prev;
        for (int iter=0; iter<MAX_ITER; ++iter) {
            std::array<long double,2> df, F, delta = {0.0,0.0};
            deriv(x_next, y_guess, df);
            computeJacobian(2, x_next, y_guess);
            for (int i=0; i<2; ++i)
                F[i] = coeff2_1[0]*y_guess[i] + coeff2_1[1]*y_prev[i] + coeff2_1[2]*y_prev2[i] - steph*df[i];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) delta[i] -= F[j]*jac[i][j];
                y_guess[i] += delta[i];
            }
            long double maxd = 0;
            for (int i=0; i<2; ++i)
                maxd = std::max(maxd, std::fabs(delta[i])/(ATOL + RTOL*std::fabs(y_guess[i]-delta[i])));
            if (maxd < 1e-1) break;
        }
        rho_chi_table[step] = y_guess[0];
        rho_cft_table[step] = y_guess[1];
        Gamma_table[step] = Gamma0 * std::pow(std::exp(x_next), Gammad);
        H_table[step]     = H(x_next, y_guess);
        y_next = y_guess;
    }

    void bdf3(long double x_next,
              const std::array<long double,2>& y_prev,
              const std::array<long double,2>& y_prev2,
              const std::array<long double,2>& y_prev3,
              std::array<long double,2>& y_next,
              int step) {
        std::array<long double,2> y_guess = y_prev;
        for (int iter=0; iter<MAX_ITER; ++iter) {
            std::array<long double,2> df, F, delta = {0.0,0.0};
            deriv(x_next, y_guess, df);
            computeJacobian(3, x_next, y_guess);
            for (int i=0; i<2; ++i)
                F[i] = coeff3_1[0]*y_guess[i] + coeff3_1[1]*y_prev[i]
                     + coeff3_1[2]*y_prev2[i] + coeff3_1[3]*y_prev3[i] - steph*df[i];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) delta[i] -= F[j]*jac[i][j];
                y_guess[i] += delta[i];
            }
            long double maxd = 0;
            for (int i=0; i<2; ++i)
                maxd = std::max(maxd, std::fabs(delta[i])/(ATOL + RTOL*std::fabs(y_guess[i])));
            if (maxd < 1.0) break;
        }
        rho_chi_table[step] = y_guess[0];
        rho_cft_table[step] = y_guess[1];
        Gamma_table[step]   = Gamma0 * std::pow(std::exp(x_next), Gammad);
        H_table[step]       = H(x_next, y_guess);
        y_next = y_guess;
    }

    void bdf4(long double x_next,
              const std::array<long double,2>& y_prev,
              const std::array<long double,2>& y_prev2,
              const std::array<long double,2>& y_prev3,
              const std::array<long double,2>& y_prev4,
              std::array<long double,2>& y_next,
              int step) {
        std::array<long double,2> y_guess = y_prev;
        for (int iter=0; iter<MAX_ITER; ++iter) {
            std::array<long double,2> df, F, delta = {0.0,0.0};
            deriv(x_next, y_guess, df);
            computeJacobian(4, x_next, y_guess);
            for (int i=0; i<2; ++i)
            {
                F[i] = coeff4_1[0]*y_guess[i] + coeff4_1[1]*y_prev[i]
                     + coeff4_1[2]*y_prev2[i] + coeff4_1[3]*y_prev3[i]
                     + coeff4_1[4]*y_prev4[i] - steph*df[i];
            }
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) delta[i] -= F[j]*jac[i][j];
                y_guess[i] += delta[i];
            }
            long double maxd = 0;
            for (int i=0; i<2; ++i)
                maxd = std::max(maxd, std::fabs(delta[i])/(ATOL + RTOL*std::fabs(y_guess[i])));
            if (maxd < 1.0) break;
        }
        rho_chi_table[step] = y_guess[0];
        rho_cft_table[step] = y_guess[1];
        Gamma_table[step]   = Gamma0 * std::pow(std::exp(x_next), Gammad);
        H_table[step]       = H(x_next, y_guess);
        y_next = y_guess;
    }
};

int main() {
    // Test inputs
    long double Gamma0_input = 1e2;   // Gyr^-1
    long double H0_input     = 70.0;  // km/s/Mpc

    Background bg;
    // Set parameters
    bg.atr        = 1e-3;
    bg.DNeff      = 0.5;
    bg.H0         = H0_input * 1e3 / c_light;
    bg.Gamma0     = Gamma0_input * Conversion * 1e3 / c_light;
    bg.Gammad     = 0.0;
    bg.Omega0_b   = 0.02238280;
    bg.Omega0_cdm = 0.3;
    bg.Omega0_chi = 0.3;
    bg.Omega0_g   = 5e-5;
    bg.Omega0_lambda = 0.67;
    bg.Omega0_ur  = 1e-6;

    if (bg.solve() != 0) {
        std::cerr << "Error solving background." << std::endl;
        return 1;
    }

    std::ofstream output("output.dat");
    output << std::setprecision(15);
    for (int i = 0; i < bg.a_size; ++i) {
        output << bg.a_table[i] << '\t'
               << bg.rho_chi_table[i] << '\t'
               << bg.rho_cft_table[i] << '\t'
               << bg.Gamma_table[i] << '\t'
               << bg.H_table[i]   << '\n';
    }
    output.close();

    std::ofstream params("params.dat");
    params << "aeq\t"    << AEQ   << '\n'
           << "Gammad\t" << bg.Gammad << '\n'
           << "Gamma0\t" << Gamma0_input << '\n'
           << "DNeff\t"  << bg.DNeff << '\n'
           << "atr\t"    << bg.atr << '\n'
           << "ainit\t"   << AINIT << '\n';
    params.close();

    std::cout << "Done." << std::endl;
    return 0;
}
