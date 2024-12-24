import numpy as np
from random_shape import tau
from parameters import F_acc, a


G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


tau_matrix = np.empty((0, (2*F_acc+1)**2))
for i in range((2*F_acc+1)**2):
    tau_row = np.array([])
    for j in range((2*F_acc+1)**2):
        vec = G_vec[i]-G_vec[j]
        tau_row = np.append(tau_row, tau((vec[0], vec[1])))
    tau_matrix = np.vstack((tau_matrix, tau_row))


tau_matrix_inverse = np.linalg.inv(tau_matrix)
np.savetxt('tau_matrix_inverse_random_shape.txt', tau_matrix_inverse)

