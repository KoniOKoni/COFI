import numpy as np
import matplotlib.pyplot as plt
import subprocess

subprocess.run(["gcc", "-o", "a.out", "boltzmann.c"])
subprocess.run(["./a.out"])

dat = open("output_backward.dat", "r")
dat_fo = open("output_forward.dat", "r")
paramdat = open("params.dat", "r")
dat.readline()
dat_fo.readline()


params = {}
a = []
rhoddm = []
rhoCFT = []
Hs = []
Gammas = []


while True:
    l = dat.readline()
    if not l:
        break
    loga, logrho_ddm, logrho_CFT, H, Gamma = list(map(float, l.strip().split()))
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
rhoddm = np.exp(np.array(rhoddm))
rhoCFT = np.exp(np.array(rhoCFT))
Hs = np.array(Hs)
Gammas = np.array(Gammas)

fig, axes = plt.subplots(2,2, figsize = (16,8), constrained_layout = True)
axes[0][0].plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0][0].plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
axes[1][0].plot(a, Gammas/Hs, 'k-', label=r"$\Gamma(a)/H(a)$")
axes[0][0].set_xscale('log')
axes[0][0].set_yscale('log')
axes[1][0].set_xscale('log')
axes[1][0].set_yscale('log')
axes[0][0].set_title(r"$\Gamma_d = {}$".format(params["Gammad"]))
axes[1][0].set_title(r"$\Gamma/H$")
axes[0][0].legend(loc=0)
axes[0][0].set_xlabel(r"$a$")
axes[1][0].set_xlabel(r"$a$")
axes[0][0].set_ylabel(r"$\rho$")
axes[1][0].set_ylabel(r"$\Gamma/H$")

a = []
rhoddm = []
rhoCFT = []
Hs = []
Gammas = []


while True:
    l = dat_fo.readline()
    if not l:
        break
    loga, logrho_ddm, logrho_CFT, H, Gamma = list(map(float, l.strip().split()))
    a.append(loga)
    rhoddm.append(logrho_ddm)
    rhoCFT.append(logrho_CFT)
    Hs.append(H)
    Gammas.append(Gamma)
    
a = np.exp(np.array(a))
rhoddm = np.exp(np.array(rhoddm))
rhoCFT = np.exp(np.array(rhoCFT))
Hs = np.array(Hs)
Gammas = np.array(Gammas)

axes[0][1].plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0][1].plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
axes[1][1].plot(a, Gammas/Hs, 'k-', label=r"$\Gamma(a)/H(a)$")
axes[0][1].set_xscale('log')
axes[0][1].set_yscale('log')
axes[1][1].set_xscale('log')
axes[1][1].set_yscale('log')
axes[0][1].set_title(r"$\Gamma_d = {}$".format(params["Gammad"]))
axes[1][1].set_title(r"$\Gamma/H$")
axes[0][1].legend(loc=0)
axes[0][1].set_xlabel(r"$a$")
axes[1][1].set_xlabel(r"$a$")
axes[0][1].set_ylabel(r"$\rho$")
axes[1][1].set_ylabel(r"$\Gamma/H$")


plt.show()
dat.close()
dat_fo.close()