from scipy.special import jv
import numpy as np
from parameters import f, a, Filling_shear, Filling_density, Background_shear, Background_density


kx = 0
ky = 0
r = (a*a*f/np.pi)**(1/2)
omega = np.empty((0, 2))


for F_acc in range(5, 31):
    omega_row = np.array([])
    G_vec = np.empty((0, 2))
    for i in range(-F_acc, F_acc + 1):
        for j in range(-F_acc, F_acc + 1):
            G_vec = np.vstack((G_vec, np.array([2 * np.pi / a * i, 2 * np.pi / a * j])))


    rho_matrix = np.empty((0, (2*F_acc+1)**2))
    for i in range((2*F_acc+1)**2):
        rho_row = np.array([])
        for j in range((2*F_acc+1)**2):
            vec = G_vec[i]-G_vec[j]
            if vec[0] == 0 and vec[1] == 0:
                rho = Filling_density*f+Background_density*(1-f)
                rho_row = np.append(rho_row, rho)
            else:
                G = (vec.dot(vec)) ** (1 / 2)
                rho = 2*(Filling_density-Background_density)*f*jv(1, G*r)/(G*r)
                rho_row = np.append(rho_row, rho)
        rho_matrix = np.vstack((rho_matrix, rho_row))

    rho_matrix_inverse = np.linalg.inv(rho_matrix)

    tau_matrix = np.empty((0, (2*F_acc+1)**2))
    for i in range((2*F_acc+1)**2):
        tau_row = np.array([])
        for j in range((2*F_acc+1)**2):
            vec = G_vec[i]-G_vec[j]
            if vec[0] == 0 and vec[1] == 0:
                tau = 1/Filling_shear*f+1/Background_shear*(1-f)
                tau_row = np.append(tau_row, tau)
            else:
                G = (vec.dot(vec)) ** (1 / 2)
                tau = 2*(1/Filling_shear-1/Background_shear)*f*jv(1, G*r)/(G*r)
                tau_row = np.append(tau_row, tau)
        tau_matrix = np.vstack((tau_matrix, tau_row))

    tau_matrix_inverse = np.linalg.inv(tau_matrix)

    tau_matrix = np.empty((0, (2*F_acc+1)**2))
    for i in range((2*F_acc+1)**2):
        tau_row = np.array([])
        for j in range((2*F_acc+1)**2):
            vec = G_vec[i]-G_vec[j]
            if vec[0] == 0 and vec[1] == 0:
                tau = Filling_shear*f+Background_shear*(1-f)
                tau_row = np.append(tau_row, tau)
            else:
                G = (vec.dot(vec)) ** (1 / 2)
                tau = 2*(Filling_shear-Background_shear)*f*jv(1, G*r)/(G*r)
                tau_row = np.append(tau_row, tau)
        tau_matrix = np.vstack((tau_matrix, tau_row))

    tau_k = np.empty((0, (2*F_acc+1)**2))
    tau_k_inverse_rule = np.empty((0, (2*F_acc+1)**2))
    for i in range((2*F_acc+1)**2):
        tau_k_row = np.array([])
        tau_k_inverse_rule_row = np.array([])
        for j in range((2*F_acc+1)**2):
            k = np.array([kx, ky])
            tau = tau_matrix[i][j]*np.dot(k+G_vec[i], k+G_vec[j])
            tau_inverse_rule = tau_matrix_inverse[i][j]*np.dot(k+G_vec[i], k+G_vec[j])
            tau_k_row = np.append(tau_k_row, tau)
            tau_k_inverse_rule_row = np.append(tau_k_inverse_rule_row, tau_inverse_rule)
        tau_k = np.vstack((tau_k, tau_k_row))
        tau_k_inverse_rule = np.vstack((tau_k_inverse_rule, tau_k_inverse_rule_row))

    multiplied = np.dot(rho_matrix_inverse, tau_k)
    multiplied_inverse_rule = np.dot(rho_matrix_inverse, tau_k_inverse_rule)
    eigenvalues = np.linalg.eigvals(multiplied)
    eigenvalues_inverse_rule = np.linalg.eigvals(multiplied_inverse_rule)
    eigenvalues = eigenvalues.real
    eigenvalues_inverse_rule = eigenvalues_inverse_rule.real
    eigenvalues.sort()
    eigenvalues_inverse_rule.sort()
    eigenvalues = eigenvalues**(1/2)
    eigenvalues_inverse_rule = eigenvalues_inverse_rule**(1/2)
    omega_row = np.append(omega_row, eigenvalues[1])
    omega_row = np.append(omega_row, eigenvalues_inverse_rule[1])

    omega = np.vstack((omega, omega_row))


    if F_acc == 10:
        print('10')
    elif F_acc == 20:
        print('20')
    elif F_acc == 25:
        print('25')
    elif F_acc == 27:
        print('27')
    elif F_acc == 29:
        print('29')

np.savetxt('deviation.txt', omega)