import numpy as np
import matplotlib.pyplot as plt
import subprocess

G = 2.75e-115
PI = np.pi

subprocess.run(["clear"])
subprocess.run(["gcc","-std=gnu11" ,"-o", "a.out", "plot.c", "-lm"])
#subprocess.run(["g++","-std=c++17","-O2", "boltzmann.cpp", "-o" ,"boltzmann"])
subprocess.run(["./a.out"])

dat = open("output.dat", "r")
#dat_fo = open("output_forward.dat", "r")
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
    loga, rho_ddm, rho_CFT, Gamma, H = list(map(float, l.strip().split()))
    a.append(loga)
    rhoddm.append(rho_ddm)
    rhoCFT.append(rho_CFT)
    Hs.append(H)
    Gammas.append(Gamma)

    
a = np.exp(np.array(a))
rhoddm = np.exp(np.array(rhoddm))
rhoCFT = np.exp(np.array(rhoCFT))
Hs = np.array(Hs)
Gammas = np.array(Gammas)

fig, axes = plt.subplots(2,1, figsize = (8,8), constrained_layout = True)
axes[0].plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0].plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
#axes[0].plot(a, rhot, 'g-', label = r"Total comoving density")
axes[1].plot(a, Gammas/Hs, 'k-', label=r"$\Gamma(a)/H(a)$")
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[0].set_title("Energy densities")
axes[1].set_title(r"$\Gamma/H$")
axes[0].legend(loc=0)
axes[0].set_xlabel(r"$a$")
axes[1].set_xlabel(r"$a$")
axes[0].set_ylabel(r"Energy densities (Mpc$^{-2}$)")
axes[1].set_ylabel(r"$\Gamma/H$")

#plt.savefig("d={}.png".format(params["Gammad"]))
plt.show()
dat.close()