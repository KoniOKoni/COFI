#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from math import pi

# import classy module
from classy import Class


common_settings = {'omega_b':0.0223828,'100*theta_s':1.0399625783900215,
                   'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842,
                   'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN'}
                   # 'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0,'non linear':'halofit'}

common_settings = {'omega_b':0.0223828,'H0':68,
                   'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842,
                   'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN'}
                   # 'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0,'non linear':'halofit'}


log10z_tr = 4
z_tr = 10**log10z_tr
a_tr = 1 / (1 + z_tr)

COFI = Class()
COFI.set(common_settings)
COFI.set({'omega_ini_dcdm':0.5, 'omega_ini_dr':10, 'Gamma_dcdm': 1e7, 'a_tr': a_tr})
COFI.compute()

print(COFI.Neff())
print(COFI.Neff_max())
print(COFI.Neff_max() - COFI.Neff())
print(COFI.omega_dcdm())
print(COFI.omega_dr())
print(COFI.omega_dcdmdr())
print(COFI.omega_dcdmdr() - COFI.omega_dr() - COFI.omega_dcdm())
print(COFI.omega_Lambda())


bak = COFI.get_background()
print(bak.keys())
z = bak['z']
rho_dr = bak['(.)rho_dr']
rho_dcdm = bak['(.)rho_dcdm']
rho_g = bak['(.)rho_g']
rho_crit = bak['(.)rho_crit']

a = 1./(1.+z)

data = open('output.dat', 'r')

a_SHK = []
rho_dcdm_SHK = []
rho_dr_SHK = []

while True:
    l = data.readline()
    if not l:
        break
    tmp1, tmp2, tmp3 = map(float, l.strip().split(','))
    a_SHK.append(np.exp(tmp1))
    rho_dcdm_SHK.append(np.exp(tmp2))
    rho_dr_SHK.append(np.exp(tmp3))
    
rho_dr_SHK = np.array(rho_dr_SHK)
rho_dcdm_SHK = np.array(rho_dcdm_SHK)
a_SHK = np.array(a_SHK)

rho_g = 5.378150968509925e-05 * np.power(a_SHK, -4.)

neff = (rho_dr_SHK) / ((7./8.) * (4./11)**(4./3.) * rho_g) + COFI.Neff()

"""
plt.plot(a, rho_dr, label='$\\rho_\mathrm{dr}^\mathrm{CLASS}$')
plt.plot(a, rho_dcdm, label='$\\rho_\mathrm{dcdm}^\mathrm{CLASS}$')
plt.plot(a_SHK, rho_dr_SHK, label='$\\rho_\mathrm{dr}^\mathrm{SHK}$')
plt.plot(a_SHK, rho_dcdm_SHK, label='$\\rho_\mathrm{dcdm}^\mathrm{SHK}$')
"""

plt.plot(a_SHK, neff)

plt.xscale('log')
plt.yscale('linear')
plt.xlabel('a')
plt.ylabel('$\\rho$')
plt.ylim(2, 5)
plt.legend()
plt.title('Background energy density')
#plt.gca().invert_xaxis()

plt.show()