import numpy as np
from random_shape import rho
from parameters import F_acc, a


G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


rho_matrix = np.empty((0, (2*F_acc+1)**2))
for i in range((2*F_acc+1)**2):
    rho_row = np.array([])
    for j in range((2*F_acc+1)**2):
        vec = G_vec[i]-G_vec[j]
        rho_row = np.append(rho_row, rho((vec[0], vec[1])))
    rho_matrix = np.vstack((rho_matrix, rho_row))


rho_matrix_inverse = np.linalg.inv(rho_matrix)
np.savetxt('rho_matrix_inverse_random_shape.txt', rho_matrix_inverse)

