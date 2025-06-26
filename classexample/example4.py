#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy import interpolate

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
print(COFI.Neff_atr())
print(COFI.Neff_max())
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
Omega_r = bak['(.)Omega_r']
rho_crit = bak['(.)rho_crit']
rho_r = Omega_r * rho_crit

neff = (rho_r - rho_g) / ((7./8.) * (4./11)**(4./3.) * rho_g)
print(neff)

neff_interp = interpolate.interp1d(z, neff, bounds_error=False, fill_value='extrapolate')
print(neff_interp(z_tr))
print(max(neff))

plt.plot(z, neff)
plt.axvline(x=z_tr, color='red', linestyle='--', label='z_tr = {:.2f}'.format(z_tr))

plt.xlim(z_tr/5, z[0])
plt.ylim(2, 5)
plt.xscale('log')
plt.yscale('linear')
plt.xlabel('z')
plt.ylabel('$\Delta N_\mathrm{eff}$')
plt.title('Neff over z')
plt.gca().invert_xaxis()

plt.show()