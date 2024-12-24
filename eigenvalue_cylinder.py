import numpy as np
from parameters import a, F_acc, f, Filling_density, Filling_shear, Background_shear, Background_density


omega = np.empty((0, 10))
rho = np.loadtxt('rho_matrix_inverse_cylinder.txt')
tau = np.loadtxt('tau_matrix_inverse_cylinder.txt')


G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


for i in range(141):
    kx = np.pi/a-np.pi/a*i/141
    ky = np.pi/a-np.pi/a*i/141
    tau_k = np.empty((0, (2 * F_acc + 1) ** 2))
    for j in range((2 * F_acc + 1) ** 2):
        tau_row = np.array([])
        for k in range((2 * F_acc + 1) ** 2):
            taug = tau[j][k]*(kx**2+(G_vec[j][0] + G_vec[k][0])*kx+(G_vec[j][0] * G_vec[k][0] + G_vec[j][1] * G_vec[k][1])+ky**2+(G_vec[j][1] + G_vec[k][1])*ky)
            tau_row = np.append(tau_row, taug)
        tau_k = np.vstack((tau_k, tau_row))
    multiplied = np.dot(rho, tau_k)
    eigenvalues = np.linalg.eigvals(multiplied)
    eigenvalues.sort()
    omega_row = np.array([])
    for i in range(10):
        omega_row = np.append(omega_row, eigenvalues[i]**(1/2))
    omega = np.vstack((omega, omega_row))

for i in range(100):
    kx = np.pi/a*i/100
    ky = 0
    tau_k = np.empty((0, (2 * F_acc + 1) ** 2))
    for j in range((2 * F_acc + 1) ** 2):
        tau_row = np.array([])
        for k in range((2 * F_acc + 1) ** 2):
            taug = tau[j][k]*(kx**2+(G_vec[j][0] + G_vec[k][0])*kx+(G_vec[j][0] * G_vec[k][0] + G_vec[j][1] * G_vec[k][1])+ky**2+(G_vec[j][1] + G_vec[k][1])*ky)
            tau_row = np.append(tau_row, taug)
        tau_k = np.vstack((tau_k, tau_row))
    multiplied = np.dot(rho, tau_k)
    eigenvalues = np.linalg.eigvals(multiplied)
    eigenvalues.sort()
    omega_row = np.array([])
    for i in range(10):
        omega_row = np.append(omega_row, eigenvalues[i]**(1/2))
    omega = np.vstack((omega, omega_row))

for i in range(100):
    kx = np.pi/a
    ky = np.pi/a*i/100
    tau_k = np.empty((0, (2 * F_acc + 1) ** 2))
    for j in range((2 * F_acc + 1) ** 2):
        tau_row = np.array([])
        for k in range((2 * F_acc + 1) ** 2):
            taug = tau[j][k]*(kx**2+(G_vec[j][0] + G_vec[k][0])*kx+(G_vec[j][0] * G_vec[k][0] + G_vec[j][1] * G_vec[k][1])+ky**2+(G_vec[j][1] + G_vec[k][1])*ky)
            tau_row = np.append(tau_row, taug)
        tau_k = np.vstack((tau_k, tau_row))
    multiplied = np.dot(rho, tau_k)
    eigenvalues = np.linalg.eigvals(multiplied)
    eigenvalues.sort()
    omega_row = np.array([])
    for i in range(10):
        omega_row = np.append(omega_row, eigenvalues[i]**(1/2))
    omega = np.vstack((omega, omega_row))

# c = ((Filling_shear*f+Background_shear*(1-f))/(Filling_density*f+Background_density*(1-f)))**(1/2)
# omega = omega*a/2/np.pi/c
omega = omega.real

np.savetxt('eigenvalue_cylinder.txt', omega)