import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter
from parameters import a, F_acc, kx, ky


eigenvector = np.loadtxt("eigenvector.txt")
G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))



fig = plt.figure()
ax = fig.add_subplot()
ax.set_aspect(1)
x = np.linspace(-a/2, a/2, 100)
y = np.linspace(-a/2, a/2, 100)
X, Y = np.meshgrid(x, y)

norm = plt.Normalize(-1.6, 1.6)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])

writer = FFMpegWriter(fps=30, metadata=dict(title='Movie', artist='Me'), bitrate=3000000)
writer.setup(fig, 'eigenvector_contour.mp4', dpi=300)

for t in range(150):
    Z = 0
    for i in range((2*F_acc+1)**2):
        coeff = eigenvector[i]
        exponential_term = np.exp(1j*(G_vec[i][0]*X+G_vec[i][1]*Y+kx*X+ky*Y-eigenvector[-1]*t/150000))
        Z += coeff*exponential_term
    Z = Z.real
    plot1 = ax.contourf(X, Y, Z, cmap='jet', norm=norm)
    writer.grab_frame()
    plot1.remove()

writer.finish()
plt.show()