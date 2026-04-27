import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energy_force.txt")

time = data[:,0]
Ekin = data[:,1]
Epot = data[:,2]
Etot = data[:,3]
T_K = data[:,4]

fig, ax = plt.subplots(ncols=2, figsize=(12, 4.5))
ax[0].plot(time, Ekin, "--", label=r"$E_\mathrm{kin}$")
ax[0].plot(time, Epot, "--", label=r"$E_\mathrm{pot}$")
ax[0].plot(time, Etot, "-.", label=r"$E_\mathrm{tot}$")
ax[0].set_xlabel("t [ps]")
ax[0].set_ylabel("E [kj/mol]") 
ax[0].legend()
ax[1].plot(time, T_K, "--", label="T")
ax[1].set_xlabel("t [ps]")
ax[1].set_ylabel("T [K]") 
plt.show()
