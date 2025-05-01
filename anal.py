import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import scipy
import scipy.integrate

H0 = 70.*1000./299792458.
cnt = 0
Gs = []
frac_chi = []
frac_CFT = []

def H(x, rho):
    rhotot = rho[0] + rho[1] + H0**2 * ((0.0223828 + 0.3)*np.exp(-3*x) + (5e-5 + 1e-6)*np.exp(-4*x) + 0.67)
    return np.sqrt(rhotot)

def deriv(x, rho, gamma0):
    return [-3*rho[0] - (gamma0/H(x, rho))*rho[0], -4*rho[1] + (gamma0/H(x, rho))*rho[0]]

start = time.time()
N = 100
tol = 1e-2
for Gamma0 in np.logspace(-1, 10, N):
    a = []
    rhoddm = []
    rhoCFT = []
    Hs = []
    Gammas = []
    subprocess.run(["clear"])
    subprocess.run(["gcc","-std=gnu11" ,"-o", "anal.out", "boltzmann.c", "-lm"])
    subprocess.run(["./anal.out", str(Gamma0), "0"])

    dat = open("output.dat", "r")
    while True:
        l = dat.readline()
        if not l:
            break
        loga, rhochi, rhocft, gamma, h = map(float, l.strip().split())
        a.append(loga)
        rhoddm.append(rhochi)
        rhoCFT.append(rhocft)
        Gammas.append(gamma)
        Hs.append(h)
    if (rhoCFT[-1] <= 0):
        sol = scipy.integrate.solve_ivp(deriv, [a[-1], np.log(1e-3)], [rhoddm[-1], rhoCFT[-1]], method='BDF', args = (Gamma0*1e3/299792458,))
        frac_chi.append(abs(sol.y[0][-1] - rhoddm[0])/rhoddm[0])
        frac_CFT.append(abs(sol.y[1][-1] - rhoCFT[0])/rhoCFT[0])
        Gs.append(Gamma0)
    
    dat.close()
end = time.time()
print("It took {} sec.".format(end - start))

plt.plot(Gs, frac_chi, 'bo', label = r"$\rho_{\chi}$ fractional difference")
plt.plot(Gs, frac_CFT, 'ro', label = r"$\rho_{\mathrm{CFT}}$ fractional difference")
plt.xscale('log')
plt.xlabel(r"$\Gamma_0$")
plt.ylabel(r"Fractional difference of $\rho(a_{tr})$")
plt.legend(loc = 0)
plt.title(r"Fractional difference of $\rho(a_{tr})$ between forward and backward integration")
plt.show()
