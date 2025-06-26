#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from math import pi

# import classy module
from classy import Class


common_settings = {'omega_b':0.0223828,'100*theta_s':1.0399625783900215,
                   'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842,
                   'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN',
                   'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0,'non linear':'halofit'}


log10z_tr = 4
z_tr = 10**log10z_tr
a_tr = 1 / (1 + z_tr)

COFI = Class()
COFI.set({'omega_ini_dcdm':0.1301075, 'Gamma_dcdm': 1e5, 'a_tr': a_tr})
COFI.compute()

print(COFI.Neff())
print(COFI.Neff_atr())

bak = COFI.get_background()
print(bak.keys())
z = bak['z']
w_dr = bak['(.)w_dr']
cs2_dr = bak['(.)cs2_dr']
w_prime_dr = bak['(.)w_prime_dr']

print(w_prime_dr)

plt.plot(z, w_dr, label='$w_\mathrm{dr}$')
plt.plot(z, cs2_dr, label='$cs^2_\mathrm{dr}$')
plt.plot(z, w_prime_dr, label='$w\'_\mathrm{dr}$')

plt.xscale('log')
plt.yscale('linear')
plt.xlabel('z')
plt.legend()
plt.title('Background equation of state for DR')
plt.gca().invert_xaxis()

plt.show()