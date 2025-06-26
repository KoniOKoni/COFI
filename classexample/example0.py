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
COFI.set({'omega_ini_dcdm':0.5, 'omega_ini_dr':10, 'Gamma_dcdm': 1e10, 'a_tr': a_tr})
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

plt.plot(a, rho_dr, label='$\\rho_\mathrm{dr}$')
plt.plot(a, rho_dcdm, label='$\\rho_\mathrm{dcdm}$')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('$\\rho$')
plt.legend()
plt.title('Background energy density')
#plt.gca().invert_xaxis()

plt.show()