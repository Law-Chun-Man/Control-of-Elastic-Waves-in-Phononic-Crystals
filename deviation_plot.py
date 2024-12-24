import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["mathtext.fontset"] = "cm"

kx = 0  # change position
ky = 0

omega = np.loadtxt('deviation.txt')
omega = omega.transpose()

fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.set_xlabel('N')
ax.set_ylabel('Reduced frequency')
plt.xlim(0, 31)
ax.set_xticks(np.arange(0, 31, 2))
ax.scatter(np.arange(5, 31, 1), omega[0], c='red', marker='x', label='Laurent Rule')
ax.scatter(np.arange(5, 31, 1), omega[1], edgecolor='blue', marker='o', facecolor='none', label='Inverse Rule')
ax.legend()
plt.title(r'''Convergence of eigenfrequency at $\overline{\Gamma}$ point with
increasing number of terms of Fourier series''')
plt.tight_layout()
plt.savefig('plot/deviation.png', dpi=600)
plt.show()
