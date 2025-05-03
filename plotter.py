import numpy as np
import matplotlib.pyplot as plt
import subprocess

G = 2.75e-115
PI = np.pi

subprocess.run(["clear"])
subprocess.run(["gcc","-std=gnu11" ,"-o", "a.out", "boltzmann.c", "-lm"])
#subprocess.run(["g++","-std=c++17","-O2", "boltzmann.cpp", "-o" ,"boltzmann"])
subprocess.run(["./a.out"])

dat = open("output.dat", "r")
#dat_fo = open("output_forward.dat", "r")
paramdat = open("params.dat", "r")
#dat_fo.readline()

params = {}
a = []
rhoddm = []
rhoCFT = []
Hs = []
Gammas = []
rhot = []


while True:
    l = dat.readline()
    if not l:
        break
    loga, logrho_ddm, logrho_CFT, Gamma, H = list(map(float, l.strip().split()))
    a.append(loga)
    rhoddm.append(logrho_ddm)
    rhoCFT.append(logrho_CFT)
    Hs.append(H)
    Gammas.append(Gamma)

while True:
    l = paramdat.readline()
    if not l:
        break
    name, value = l.strip().split()
    params[name] = float(value)

    
a = np.exp(np.array(a))
rhoddm = (np.array(rhoddm))
rhoCFT = (np.array(rhoCFT))
Hs = np.array(Hs)
Gammas = np.array(Gammas)

rhoddminit = 0.3*(70*1000/299792458)**2 * np.power(a[0], -3) * params["fddm"]



print("rhoddm = {}, rhoCFT = {}".format(rhoddm[-1], rhoCFT[-1]))

fig, axes = plt.subplots(2,1, figsize = (16,8), constrained_layout = True)
axes[0].plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0].plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
axes[0].plot(a, rhoddminit*np.power(a[0]/a, 3), 'k--', label = r"$\rho_{\mathrm{ddm}}$, no decay")
#axes[0].plot(a, rhot, 'g-', label = r"Total comoving density")
axes[1].plot(a, Gammas/Hs, 'k-', label=r"$\Gamma(a)/H(a)$")
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[0].set_title(r"$\Gamma_d = {}$, $\tau = {}$ years".format(params["Gammad"], params["tau"]))
axes[1].set_title(r"$\Gamma/H$")
axes[0].legend(loc=0)
axes[0].set_xlabel(r"$a$")
axes[1].set_xlabel(r"$a$")
axes[0].set_ylabel(r"Energy densities (Mpc$^{-2}$)")
axes[1].set_ylabel(r"$\Gamma/H$")
axes[0].axvline(x=params['aeq'], color = "black", linestyle = 'dashed')
axes[0].axvline(x=9e-4, color = "green", linestyle = 'dashed')
axes[1].axvline(x=params['aeq'], color = "black", linestyle = 'dashed')
axes[1].axvline(x=9e-4, color = "green", linestyle = 'dashed')

for i in range(2):
    axes[i].set_xlim([params["ainit"], 1])


#plt.savefig("d={}.png".format(params["Gammad"]))
plt.show()
dat.close()