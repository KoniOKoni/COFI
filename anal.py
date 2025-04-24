import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time

X = []
Y = []

start = time.time()
subprocess.run(["clear"])
subprocess.run(["gcc","-std=gnu11" ,"-o", "anal.out", "anal.c", "-lm"])
subprocess.run(["./anal.out"])

dat = open("scan.dat", "r")
while True:
    l = dat.readline()
    if not l:
        break
    Gamma0, Gammad = map(float, l.strip().split())
    X.append(Gamma0)
    Y.append(Gammad)
end = time.time()
print("Takes {} s".format(end - start))
dat.close()

plt.scatter(X,Y)
plt.xscale("log")
plt.xlabel(r"$\Gamma_0$ (km/s/Mpc)")
plt.ylabel(r"$\Gamma_d$")
plt.show()
