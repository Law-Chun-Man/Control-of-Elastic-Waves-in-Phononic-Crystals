import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from matplotlib.animation import FFMpegWriter
from parameters import a, F_acc, kx, ky


eigenvector = np.loadtxt("eigenvector.txt")
G_vec = np.empty((0, 2))
for i in range(-F_acc, F_acc+1):
    for j in range(-F_acc, F_acc+1):
        G_vec = np.vstack((G_vec, np.array([2*np.pi/a*i, 2*np.pi/a*j])))



fig = plt.figure(figsize=(6.4, 3.6))
fig.set_facecolor('black')
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('black')
ax.margins(x=0, y=0, z=0)
ax.set_proj_type('ortho')
ax.set_zlim(-1.6, 1.6)
ax.tick_params(axis='z', colors='none')
ax.tick_params(axis='y', colors='none')
ax.tick_params(axis='z', colors='white')
ax.set_zlabel("displacement", c="white")
ax.set_box_aspect((1, 1, 0.5))
ax.set_position([0, 0, 1, 1])
ax.grid(False)
plt.tight_layout(pad=0)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
x = np.linspace(-a/2, a/2, 200)
y = np.linspace(-a/2, a/2, 200)
X, Y = np.meshgrid(x, y)

norm = plt.Normalize(-1.6, 1.6)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])


image = Image.open("plot/primitive.png")
pixels = image.load()
points = np.empty((0, 3))
for i in range(192):
    for j in range(192):
        if pixels[i, j] != (0, 0, 0, 255):
            points = np.vstack((points, np.array([-a*i/191+a/2, a*j/191-a/2, -1.5])))

points = np.transpose(points)
words = ax.scatter(points[0], points[1], points[2], c='cyan', s=0.05)


writer = FFMpegWriter(fps=30, metadata=dict(title='Movie', artist='Me'), bitrate=3000000)
writer.setup(fig, 'plot/eigenvector_plot.mp4', dpi=300)

for t in range(200):
    Z = 0
    for i in range((2*F_acc+1)**2):
        coeff = eigenvector[i]
        exponential_term = np.exp(1j*(G_vec[i][0]*X+G_vec[i][1]*Y+kx*X+ky*Y-eigenvector[-1]*t/150000))
        Z += coeff*exponential_term
    Z = Z.real
    plot1 = ax.plot_surface(X, Y, Z, cmap='jet', norm=norm)
    ax.view_init(elev=45-45*(t/199)**(1/2), azim=91+40*t/199)
    words.set_alpha(1-(t/199)**(1/2))
    writer.grab_frame()
    plot1.remove()

for t in range(200):
    Z = 0
    for i in range((2*F_acc+1)**2):
        coeff = eigenvector[i]
        exponential_term = np.exp(1j*(G_vec[i][0]*X+G_vec[i][1]*Y+kx*X+ky*Y-eigenvector[-1]*t/150000-eigenvector[-1]*200/150000))
        Z += coeff*exponential_term
    Z = Z.real
    plot1 = ax.plot_surface(X, Y, Z, cmap='jet', norm=norm)
    ax.view_init(elev=0, azim=131+40*t/199)
    writer.grab_frame()
    plot1.remove()

writer.finish()
plt.show()
