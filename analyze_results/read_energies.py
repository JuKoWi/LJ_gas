import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energy_force.txt")

time = data[:,0]
Ekin = data[:,1]
Epot = data[:,2]
Etot = data[:,3]

fig, ax = plt.subplots()
ax.plot(time, Ekin, label=r"$E_\mathrm{kin}$")
ax.plot(time, Epot, label=r"$E_\mathrm{pot}$")
ax.plot(time, Etot, label=r"$E_\mathrm{tot}$")
ax.legend()
plt.show()
