import numpy as np
import matplotlib.pyplot as plt
import subprocess

def cftscale(a, atr):
    if a > np.log(atr):
        return np.exp(3*a)
    else:
        return np.exp(4*a)
    
np.vectorize(cftscale)

G = 2.75e-115
PI = np.pi

subprocess.run(["clear"])
subprocess.run(["gcc", "-o", "a.out", "boltzmann.c", "-lm"])
subprocess.run(["./a.out"])

dat = open("output.dat", "r")
#dat_fo = open("output_forward.dat", "r")
paramdat = open("params.dat", "r")
dat.readline()
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
rhoddm = np.exp(np.array(rhoddm))
rhoCFT = np.exp(np.array(rhoCFT))
Hs = np.array(Hs)
Gammas = np.array(Gammas)
rhot = np.exp(np.array(rhot))

print("rhoddm = {}, rhoCFT = {}".format(rhoddm[-1], rhoCFT[-1]))

fig, axes = plt.subplots(2,1, figsize = (16,8), constrained_layout = True)
axes[0].plot(a, rhoddm, 'r-', label = r"$\rho_{\mathrm{ddm}}$")
axes[0].plot(a, rhoCFT, 'b-', label = r"$\rho_{\mathrm{CFT}}$")
#axes[0].plot(a, rhot, 'g-', label = r"$\rho_{\mathrm{tot}}$")
axes[1].plot(a, Gammas/Hs, 'k-', label=r"$\Gamma(a)/H(a)$")
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[0].set_title(r"$\Gamma_d = {}$, $\Gamma_0 = {}$ km/s/Mpc".format(params["Gammad"], float(params["Gamma0"])))
axes[1].set_title(r"$\Gamma/H$")
axes[0].legend(loc=0)
axes[0].set_xlabel(r"$a$")
axes[1].set_xlabel(r"$a$")
axes[0].set_ylabel(r"$\frac{8\pi G}{3}\rho$ (Mpc$^{-2}$)")
axes[1].set_ylabel(r"$\Gamma/H$")
axes[0].axvline(x=params['aeq'], color = "black", linestyle = 'dashed')
axes[0].axvline(x=params['atr']*0.99, color = "green", linestyle = 'dashed')
axes[1].axvline(x=params['aeq'], color = "black", linestyle = 'dashed')
axes[1].axvline(x=params['atr']*0.99, color = "green", linestyle = 'dashed')

for i in range(2):
    axes[i].set_xlim([params["ainit"], 1])


#plt.savefig("d={}.png".format(params["Gammad"]))
plt.show()
dat.close()