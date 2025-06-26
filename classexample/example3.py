#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pandas as pd

# import classy module
from classy import Class

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

COFI = Class()
COFI.set(common_settings)
COFI.set({'omega_ini_dcdm':0.1301075, 'Gamma_dcdm': 1e-1, 'a_tr': a_tr})
COFI.compute()

# CMB
kk = np.logspace(-4,np.log10(3),1000)
colors = ['tab:red', 'tab:blue', 'tab:green']
labels = ['$\mathcal{D}$', '$\mathcal{DH}$', '$\mathcal{DHF}$']

set_linewidth = 2
tick_size = 14

# plot TT    
plt.figure(1)
plt.xscale('linear');plt.yscale('linear')
plt.xlabel(r'$\ell$', fontdict=font)
plt.ylabel('$\Delta C_\ell^{\mathrm{TT}} / C^{\mathrm{TT}}_{\ell\Lambda\mathrm{CDM}}$', fontdict=font)
plt.axhline(0, color='tab:gray', linewidth=1)

# plot EE    
plt.figure(2)
plt.xscale('linear');plt.yscale('linear')
plt.xlabel(r'$\ell$', fontdict=font)
plt.ylabel('$\Delta C_\ell^{\mathrm{EE}} / C^{\mathrm{EE}}_{\ell\Lambda\mathrm{CDM}}$', fontdict=font)
plt.axhline(0, color='tab:gray', linewidth=1)


cls_lcdm = LambdaCDM.lensed_cl(2500)
ll = cls_lcdm['ell'][2:]
cls = COFI.lensed_cl(2500)

# plot TT    
plt.figure(1)
cl0 = cls_lcdm['tt'][2:]
cl = cls['tt'][2:]
plt.plot(ll, (cl-cl0)/cl0, linewidth=set_linewidth)    

# plot EE    
plt.figure(2)
cl0 = cls_lcdm['ee'][2:]
cl = cls['ee'][2:]
plt.plot(ll, (cl-cl0)/cl0, linewidth=set_linewidth)    


plt.figure(1)
legend = plt.legend(loc='best', fontsize=13)
for words in legend.get_texts():
    words.set_fontname(fnm)

plt.xticks(**fontdict)
plt.yticks(**fontdict)
plt.xlim(20,2500); plt.ylim(-1, 1); 

plt.xlabel(r'$\ell$', fontdict=font)
plt.ylabel('$\Delta C_\ell^{\mathrm{TT}} / C^{\mathrm{TT}}_{\ell\Lambda\mathrm{CDM}}$', fontdict=font)
plt.title(r'TT residuals for best-fits', **fontdict)

plt.savefig('./plots/bf_TT.pdf', bbox_inches="tight")
plt.close()

plt.figure(2)
legend = plt.legend(loc='best', fontsize=13)
for words in legend.get_texts():
    words.set_fontname(fnm)


plt.xticks(**fontdict)
plt.yticks(**fontdict)
plt.xlim(20,2500); plt.ylim(-1, 2); 

plt.xlabel(r'$\ell$', fontdict=font)
plt.ylabel('$\Delta C_\ell^{\mathrm{EE}} / C^{\mathrm{EE}}_{\ell\Lambda\mathrm{CDM}}$', fontdict=font)
plt.title(r'EE residuals for best-fits', **fontdict)

plt.savefig('./plots/bf_EE.pdf', bbox_inches="tight")
plt.close()
