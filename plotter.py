import numpy as np
import matplotlib.pyplot as plt
import subprocess


subprocess.run(["clear"])
subprocess.run(["rm", "-rf", "plots/*.png"])
subprocess.run(["gcc","-std=gnu11" ,"-o", "a.out", "plot.c", "-lm"])
subprocess.run(["./a.out"])

omega_table = np.logspace(-2, np.log10(0.5), num=50 ,endpoint=True)

for i in range(50):
    print("Plotting {}-th plot.".format(i+1))
    subprocess.run(["./a.out", "{}".format(omega_table[i])])
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
    plt.figure()
    plt.plot(a, rhoddm, label = r"$\rho_{dcdm}$")
    plt.plot(a, rhoCFT, label = r"$\rho_{dr}$")
    plt.legend(loc = 0)
    plt.title(r'BED, $\Omega_{{dcdm, ini}}$ = {}, $\Gamma = 10^{{6}}$ km/s/Mpc'.format(omega_table[i]))
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("plots/{}.png".format(i))
    plt.close()
    dat.close()
