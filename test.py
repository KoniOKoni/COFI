import numpy as np
import matplotlib.pyplot as plt

out = open("output.dat", 'r')
tout = open('test.dat', 'r')

a = []
rho = []
cft = []
at = []
rhot = []
cftt = []
h = []
g = []

while True:
    l = out.readline()
    if not l:
        break
    loga, logddm, logCFT, H, gamma = map(float, l.strip().split())
    a.append(loga)
    rho.append(logddm)
    cft.append(logCFT)
    h.append(H)
    g.append(gamma)

while True:
    l = tout.readline()
    if not l:
        break
    loga, logddm, logCFT= map(float, l.strip().split())
    at.append(loga)
    rhot.append(logddm)
    cftt.append(logCFT)

a = np.exp(np.array(a))
rho = np.exp(np.array(rho))
cft = np.exp(np.array(cft))
at = np.exp(np.array(at))
rhot = np.exp(np.array(rhot))
cftt = np.exp(np.array(cftt))
g = np.array(g)
h = np.array(h)

fig, axes = plt.subplots(2,1, figsize = (16,8))

axes[0].set_title(r"$\Gamma_d = -3.5$")
axes[0].plot(a, rho, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0].plot(a, cft, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
axes[0].plot(at, rhot, 'ko', markersize = 2)
axes[0].plot(at, cftt, 'go', markersize = 2)
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlabel(r"$a$")
axes[0].set_ylabel(r"$\rho$")
axes[0].vlines(2.94e-4, 0, max(rho[-1], cft[-1]), colors='black', linestyles='dashed')
axes[0].legend(loc=0)

axes[1].plot(a, g/h, 'k-', label = r"$\Gamma/H$")
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel(r"$a$")
axes[1].set_ylabel(r"$\Gamma/H$")
axes[1].legend(loc=0)
axes[1].vlines(2.94e-4, 0, max(g/h), colors='black', linestyles='dashed')
plt.show()

out.close()
tout.close()