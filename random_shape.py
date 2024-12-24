from PIL import Image
import numpy as np
from parameters import acc, a, Filling_shear, Filling_density, Background_shear, Background_density
import matplotlib.pyplot as plt

image = Image.open("plot/heart.png")
pixels = image.load()


def shear(r):
    x, y = r
    inside_cylinder = 1

    if pixels[x, y] == (255, 255, 255, 255):
        inside_cylinder = 0  # check if the input point is in the cylinder

    if inside_cylinder == 1:
        return 1/Filling_shear
    if inside_cylinder == 0:
        return 1/Background_shear


def tau(G):
    x, y = G
    sum = 0
    for i in range(acc+1):
        for j in range(acc+1):
            i_1 = a*i/acc-a/2
            j_1 = a*j/acc-a/2
            c = shear((i, j))*np.exp(-1j*(i_1*x+j_1*y))*(a/acc)**2
            sum += c  # sum over the area of a unit cell using numerical integration
    return (sum/a**2)


def density(r):
    x, y = r
    inside_cylinder = 1

    if pixels[x, y] == (255, 255, 255, 255):
        inside_cylinder = 0  # check if the input point is in the cylinder

    if inside_cylinder == 1:
        return Filling_density
    if inside_cylinder == 0:
        return Background_density


def rho(G):
    x, y = G
    sum = 0
    for i in range(acc+1):
        for j in range(acc+1):
            i_1 = a*i/acc-a/2
            j_1 = a*j/acc-a/2
            c = density((i, j))*np.exp(-1j*(i_1*x+j_1*y))*(a/acc)**2
            sum += c  # sum over the area of a unit cell using numerical integration
    return (sum/a**2)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X, Y, Z = np.array([]), np.array([]), np.array([])
# filling_count = 0
# background_count = 0
# for i in range(acc+1):
#     for j in range(acc+1):
#         X = np.append(X, i)
#         Y = np.append(Y, j)
#         Z = np.append(Z, shear((i, j)))
#         if shear((i, j)) == 1/Filling_shear:
#             filling_count += 1
#         elif shear((i, j)) == 1/Background_shear:
#             background_count += 1
#
# ax.scatter(X, Y, Z, s=1)
# plt.show()
# print(filling_count/(filling_count+background_count))
# print(filling_count+background_count)
