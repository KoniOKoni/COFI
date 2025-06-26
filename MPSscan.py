#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pandas as pd
import subprocess

# import classy module
from classy import Class

subprocess.run(["rm", "-rf", "plots"])
subprocess.run(["mkdir", "plots"])


plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'#['dejavusans', 'dejavuserif', 'cm', 'stix', 'stixsans', 'custom']
font = {'family': 'Times New Roman',
        #'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
fontdict = {'fontsize':14.,
            'fontname':'Times New Roman'}
fnm = 'Times New Roman'

# Lambda CDM
LambdaCDM = {}
var_data = ["./D1_LCDM.txt"]
for idx, val in enumerate(var_data):
    data = pd.read_csv(val, sep="\t", header=None, dtype=str).to_numpy()
    data = data[0]
    print(data)
    common_settings = {'omega_b':data[0], '100*theta_s':data[2],
                        'ln10^{10}A_s':data[3],'n_s':data[4],'tau_reio':data[5],
                        'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN',
                    'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_h/Mpc':3.0,'non linear':'halofit'}

    LambdaCDM = Class()
    LambdaCDM.set(common_settings)
    LambdaCDM.compute()


# CDR
common_settings = {'omega_b':0.0223828,'100*theta_s':1.0399625783900215,
                   'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842,
                   'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN',
                   'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0,'non linear':'halofit'}


log10z_tr = 4
z_tr = 10**log10z_tr
a_tr = 1 / (1 + z_tr)

dcdm_table = np.logspace(-2, 0, 10)
dr_table = np.logspace(-2, 1, 10)

# CMB
kk = np.logspace(-4,np.log10(3),1000)
colors = ['tab:red', 'tab:blue', 'tab:green']
labels = ['$\mathcal{D}$', '$\mathcal{DH}$', '$\mathcal{DHF}$']
pkM_LCDM = []

for k in kk:
    pkM_LCDM.append(LambdaCDM.pk(k, 0.))
    
set_linewidth = 2
tick_size = 14
cnt = 1

for i in range(10):
    for j in range(10):
        pkM = []
        COFI = Class()
        COFI.set(common_settings)
        COFI.set({'omega_ini_dcdm':dcdm_table[i], 'omega_ini_dr':dr_table[j], 'Gamma_dcdm': 1e6, 'a_tr': a_tr})
        COFI.compute()
        
        for k in kk:
            pkM.append(COFI.pk(k,0.))

        # plot TT    
        plt.figure(1)
        plt.xscale('log');plt.yscale('linear')
        plt.xlabel(r'$k$', fontdict=font)
        plt.ylabel('$P(k)$', fontdict=font)
        plt.axhline(0, color='tab:gray', linewidth=1)

        # plot TT    
        plt.figure(1)
        cl0 = np.array(pkM_LCDM)
        cl = np.array(pkM)
        plt.plot(kk, (cl - cl0)/cl0, linewidth=set_linewidth)      

        plt.figure(1)
        legend = plt.legend(loc='best', fontsize=13)
        for words in legend.get_texts():
            words.set_fontname(fnm)

        plt.xticks(**fontdict)
        plt.yticks(**fontdict)
        plt.xlim(1e-4,3); plt.ylim(-2, 5); 

        plt.xlabel(r'$k$', fontdict=font)
        plt.ylabel(r'$\Delta P(k)/P(k)_{\Lambda \mathrm{CDM}}$', fontdict=font)
        plt.title(r'MPS residuals, $\Omega_{{dcdm, ini}}$ = {:.2f}, $\Omega_{{dr, ini}}$ = {:.2f} $\Gamma = 10^{{6}}$ km/s/Mpc'.format(dcdm_table[i], dcdm_table[j]), **fontdict)

        plt.savefig('plots/MPS_{}.pdf'.format(cnt), bbox_inches="tight")
        plt.close()

        cnt += 1
