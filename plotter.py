import numpy as np
import matplotlib.pyplot as plt

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
    
a = np.array(a)
rhoddm = np.array(rhoddm)
rhoCFT = np.array(rhoCFT)

plt.plot(np.exp(a), np.exp(rhoddm), 'r-', label = r"$\rho_{\mathrm{ddm}}$")
plt.plot(np.exp(a), np.exp(rhoCFT), 'b-', label = r"$\rho_{\mathrm{CFT}}$")
plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 0)
plt.xlabel(r"$a$")
plt.ylabel(r"$\rho$ (eV$^4$)")
plt.show()
dat.close()