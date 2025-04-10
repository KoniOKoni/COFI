import numpy as np
import matplotlib.pyplot as plt
import subprocess

subprocess.run(["gcc", "-o", "a.out", "boltzmann.c"])
subprocess.run(["./a.out"])

dat = open("output.dat", "r")
dat.readline()

a = []
rhoddm = []
rhoCFT = []


while True:
    l = dat.readline()
    if not l:
        break
    loga, logrho_ddm, logrho_CFT = list(map(float, l.strip().split()))
    a.append(loga)
    rhoddm.append(logrho_ddm)
    rhoCFT.append(logrho_CFT)
    
a = np.exp(np.array(a))
rhoddm = np.exp(np.array(rhoddm))
rhoCFT = np.exp(np.array(rhoCFT))

plt.plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
plt.plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 0)
plt.vlines(3e-4, 0, max(rhoddm[-1], rhoCFT[-1]), colors = 'black', linestyles='dashed')
plt.xlabel(r"$a$")
plt.ylabel(r"$\rho$ (eV$^4$)")
plt.show()
dat.close()