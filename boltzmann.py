import numpy as np
import matplotlib.pyplot as plt
import scipy
from math import *

#To follow CLASS convention, every quantitis are expressed in natural units
#in which everything is in powers of Mpc

#Conversion factors
Mpc_over_m =  3.085677581282e22
Gyr_over_Mpc = 3.06601394e2
c = 299792458 #Speed of light in m/s

#Global inputs. Omegas are today's values.
ainit = 1e-7
aeq = 3e-4
atr = 1e-3
DNeff = 1.0
Gamma0 = 1e3 #Gyr^-1
Gammad = -2
Omega_ddm = 1e-10
Omega_b = 0.049389
Omega_cdm = 0.265028
Omega_Lambda = 0.679
Omega_photon = 8.4e-5
Omega_neutrino = 1e-6
H0 = 67.32 #km * s^-1 * Mpc^-1
PI = np.pi

#Convert into powers of Mpc
H0 = 1000*H0/c #Mpc^-1
Gamma0 = Gamma0*Gyr_over_Mpc #Mpc^-1
G = 2.75e-115 #Mpc^2
rhocrit = 2.19e106 #Mpc^-4

rho_CFT_atr = H0**2 * Omega_photon * pow(atr, -4) * (7.0/8.0) * pow(4.0/11.0, 4.0/3.0) * DNeff
rho_ddm_atr = H0**2 * Omega_ddm * pow(atr, -3)



#Conformal Hubble constant
def H(x, y):
    rhotot = y[0] + y[1] + rhocrit*((Omega_neutrino + Omega_photon)*exp(-4*x) + (Omega_b + Omega_cdm)*exp(-3*x) + Omega_Lambda)
    if rhotot < 0:
        print("At {}, rhotot is negative = {}".format(x, rhotot))
        return 1e-100
    return sqrt(rhotot)

#y = (rho_ddm, rho_CFT)
#x = log(a)
def BoltzmannEQ(x, y):
    drhoddmdx = -3*y[0] - Gamma0/H(x,y) * exp((Gammad*x)) * y[0]
    drhoCFTdx = -4*y[1] + Gamma0/H(x,y) * exp((Gammad*x)) * y[0]
    return [drhoddmdx, drhoCFTdx]

#Solve and plot
sol = scipy.integrate.solve_ivp(BoltzmannEQ, [log(atr), log(ainit)], [rho_ddm_atr, rho_CFT_atr], dense_output = True, method = 'BDF', rtol = 1e-8, atol = 1e-12)
a = np.linspace(log(atr), log(ainit), 1000)
sol = sol.sol(a)
plt.plot(a, sol[0], 'r-', label = r"$\rho_{\mathrm{ddm}}$")
plt.plot(a, sol[1], 'b-', label = r"$\rho_{\mathrm{CFT}}$")
if aeq < atr:
    plt.vlines(log(aeq), 0, max(sol[0][-1], sol[1][-1]), colors = 'black', linestyles='dashed')
    
plt.yscale('log')
plt.xlabel(r"$\log(a)$")
plt.ylabel(r"$\frac{8\pi G}{3}\rho(a)$ (Mpc$^{-2}$)")
plt.legend(loc = 0)
plt.show()

