import numpy as np
from parameters import a, F_acc, kx, ky


G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


rho_matrix_inverse = np.loadtxt("rho_matrix_inverse_cylinder.txt", dtype=complex)
tau_matrix_inverse = np.loadtxt("tau_matrix_inverse_cylinder.txt", dtype=complex)


tau_k = np.empty((0, (2 * F_acc + 1) ** 2))
for i in range((2 * F_acc + 1) ** 2):
    tau_k_row = np.array([])
    for j in range((2 * F_acc + 1) ** 2):
        k = np.array([kx, ky])
        tau = tau_matrix_inverse[i][j] * np.dot(k + G_vec[i], k + G_vec[j])
        tau_k_row = np.append(tau_k_row, tau)
    tau_k = np.vstack((tau_k, tau_k_row))


multiplied = np.dot(rho_matrix_inverse, tau_k)
eigenvalues, eigenvectors = np.linalg.eig(multiplied)
eigenvectors = np.transpose(eigenvectors)
index = np.argsort(eigenvalues)
eigenvector = eigenvectors[index[0]]
eigenvector = np.append(eigenvector, eigenvalues[index[0]]**(1/2))
eigenvector = eigenvector.real
np.savetxt("eigenvector.txt", eigenvector)