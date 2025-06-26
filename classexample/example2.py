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
var_k = 0.5

COFI = Class()
COFI.set(common_settings)
COFI.set({'omega_ini_dcdm':0.13, 'Gamma_dcdm': 1e4, 'a_tr': a_tr})
COFI.set({'k_output_values':var_k})
COFI.compute()

all_k = COFI.get_perturbations()  # this potentially constains scalars/tensors and all k values
print(all_k['scalar'][0].keys())
one_k = all_k['scalar'][0]
delta_dm = one_k['delta_dcdm']
delta_g = one_k['delta_g']
delta_dr = one_k['delta_dr']
tau = one_k['tau [Mpc]']
a = one_k['a']

plt.plot(a, np.abs(delta_dm), label='$\delta_\mathrm{dm}$')
plt.plot(a, np.abs(delta_g), label='$\delta_g$')
plt.plot(a, np.abs(delta_dr), label='$\delta_\mathrm{dr}$')

plt.yscale('log')
plt.xscale('log')
plt.xlabel('$a$')
plt.ylabel('$\delta$')
plt.legend(loc='best')

plt.show()