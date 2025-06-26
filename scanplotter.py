import numpy as np
import matplotlib.pyplot as plt
import subprocess
from classy import Class


subprocess.run(["clear"])
subprocess.run(["rm", "-f", "plots/*.pdf"])
subprocess.run(["gcc","-std=gnu11" ,"-o", "a.out", "plot.c", "-lm"])
subprocess.run(["./a.out"])

dcdm_table = np.logspace(-2, 0, 10)
dr_table = np.logspace(-2, 1, 10)

common_settings = {'omega_b':0.0223828,'H0':68,
                   'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842,
                   'N_ur':2.046,'N_ncdm':1,'m_ncdm':0.06,'T_ncdm':0.7137658555036082,'YHe':'BBN'}
                   # 'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0,'non linear':'halofit'}


log10z_tr = 4
z_tr = 10**log10z_tr
a_tr = 1 / (1 + z_tr)
cnt = 1

for i in range(10):
    for j in range(10):
        COFI = Class()
        COFI.set(common_settings)
        COFI.set({'omega_ini_dcdm':dcdm_table[i], 'omega_ini_dr':dr_table[j], 'Gamma_dcdm': 1e6, 'a_tr': a_tr})
        COFI.compute()
        print("Plotting {}-th plot.".format(cnt))
        subprocess.run(["./a.out", "{}".format(dcdm_table[i]), "{}".format(dr_table[j])])
        dat = open("output.dat", "r")
        a = []
        rhoddm = []
        rhoCFT = []
        while True:
            l = dat.readline()
            if not l:
                break
            loga, rho_ddm, rho_CFT = list(map(float, l.strip().split(',')))
            a.append(loga)
            rhoddm.append(rho_ddm)
            rhoCFT.append(rho_CFT)
    
        a = np.exp(np.array(a))
        rhoddm = np.exp(np.array(rhoddm))
        rhoCFT = np.exp(np.array(rhoCFT))
        bak = COFI.get_background()
        z = bak['z']
        rho_dr = bak['(.)rho_dr']
        rho_dcdm = bak['(.)rho_dcdm']
        rho_g = bak['(.)rho_g']
        rho_crit = bak['(.)rho_crit']

        a_class = 1./(1.+z)
        plt.figure()
        plt.plot(a, rhoddm, label = r"$\rho_{dcdm}^{\mathrm{SHK}}$")
        plt.plot(a, rhoCFT, label = r"$\rho_{dr}^{\mathrm{SHK}}$")
        plt.plot(a_class, rho_dr, label='$\\rho_\mathrm{dr}^{\mathrm{CLASS}}$')
        plt.plot(a_class, rho_dcdm, label='$\\rho_\mathrm{dcdm}^{\mathrm{CLASS}}$')
        plt.legend(loc = 0)
        plt.title(r'BED, $\Omega_{{dcdm, ini}}$ = {}, $\Omega_{{dr, ini}}$ = {} $\Gamma = 10^{{6}}$ km/s/Mpc'.format(dcdm_table[i], dr_table[j]))
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig("plots/{}.png".format(cnt))
        plt.close()
        dat.close()
        cnt += 1