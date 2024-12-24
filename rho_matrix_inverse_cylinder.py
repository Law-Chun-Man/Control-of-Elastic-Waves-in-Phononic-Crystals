from scipy.special import jv
import numpy as np
from parameters import a, F_acc, f, Filling_density, Background_density


G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


r = (a*a*f/np.pi)**(1/2)
rho_matrix = np.zeros(((2*F_acc+1)**2, (2*F_acc+1)**2))
for i in range((2*F_acc+1)**2):
    for j in range((2*F_acc+1)**2):
        vec = G_vec[i]-G_vec[j]
        if vec[0] == 0 and vec[1] == 0:
            rho = Filling_density*f+Background_density*(1-f)
            rho_matrix[i][j] = rho
        else:
            G = vec
            G = (G.dot(G)) ** (1 / 2)
            rho = 2*(Filling_density-Background_density)*f*jv(1, G*r)/(G*r)
            rho_matrix[i][j] = rho

rho_matrix_inverse = np.linalg.inv(rho_matrix)
np.savetxt('rho_matrix_inverse_cylinder.txt', rho_matrix_inverse)
