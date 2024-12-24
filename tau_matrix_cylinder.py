from scipy.special import jv
import numpy as np
from parameters import a, F_acc, f, Filling_shear, Background_shear

G_vec = np.empty((0, 2))

for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))


r = (a*a*f/np.pi)**(1/2)
tau_matrix = np.zeros(((2*F_acc+1)**2, (2*F_acc+1)**2))
for i in range((2*F_acc+1)**2):
    for j in range((2*F_acc+1)**2):
        vec = G_vec[i]-G_vec[j]
        if vec[0] == 0 and vec[1] == 0:
            tau = Filling_shear*f+Background_shear*(1-f)
            tau_matrix[i][j] = tau
        else:
            G = vec
            G = (G.dot(G)) ** (1 / 2)
            tau = 2*(Filling_shear-Background_shear)*f*jv(1, G*r)/(G*r)
            tau_matrix[i][j] = tau


np.savetxt('tau_matrix.txt', tau_matrix)