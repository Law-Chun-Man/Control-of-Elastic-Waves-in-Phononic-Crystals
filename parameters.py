import numpy as np

F_acc = 10 # Fourier series half length
acc = 100  # define the accuracy of the numerical integration
f = 0.6  # fraction of cylinder in the cell
a = 1  # side length of the lattice
Filling_shear = 8.847e10  # rho*c_t^2 of filling material in SI unit
Background_shear = 1.615e9  # rho*c_t^2 of background material in SI unit
Filling_density = 1750  # density of filling material in SI unit
Background_density = 1200  # density of background material in SI unit
kx = np.pi/a
ky = np.pi/a