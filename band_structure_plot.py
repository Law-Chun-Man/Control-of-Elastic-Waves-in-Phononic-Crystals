import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from parameters import a, f
plt.rcParams["mathtext.fontset"] = "cm"


omega = np.loadtxt('eigenvalue_cylinder.txt')
omegat = np.transpose(omega)


fig, ax = plt.subplots(figsize=(6, 7))

ax.set_xlabel('Reduced wave vector')
ax.set_ylabel('Reduced frequency')
plt.xlim(-2**(1/2)*np.pi/a, 2*np.pi/a)
ax.set_xticks([-2**(1/2)*np.pi/a, 0, np.pi/a, 2*np.pi/a])
ax.set_xticklabels([r'$\overline{\mathrm{M}}$', r'$\overline{\Gamma}$', r'$\overline{\mathrm{X}}$', r'$\overline{\mathrm{M}}$'], fontsize=16)
plt.ylim(0, 2)
ax.set_yticks([0, 0.5, 1, 1.5, 2])
ax.plot([0, 0], [0, 2], c='black')
ax.plot([np.pi/a, np.pi/a], [0, 2], c='black')
plt.title('Band structure')



for i in range(10):
    ax.plot(np.linspace(-2**(1/2)*np.pi/a, 2*np.pi/a, len(omega)), omegat[i], c='black')


# band gap
gap = np.array([])
for i in range(len(omegat)-1):
    max = np.max(omegat[i])
    min = np.min(omegat[i+1])
    if max < min:
        gap = np.append(gap, max)
        gap = np.append(gap, min)
for i in range(0, len(gap), 2):
    lower = gap[i]
    upper = gap[i+1]
    ax.fill_between([-2 ** (1 / 2) * np.pi / a, 2 * np.pi / a], lower, upper, color='none', hatch='///',
                    edgecolor='black')


r = (a*a*f/np.pi)**(1/2)
rectangle = patches.Rectangle((2*np.pi/a-2.25, 0), 2.25, 0.36, edgecolor='black', facecolor='#7D2882')
ax.add_patch(rectangle)
ellipse = patches.Ellipse((2*np.pi/a-2.25/2, 0.36/2), r/0.5*2.25, r/0.5*0.36, edgecolor='none', facecolor='#F0AA23')
ax.add_patch(ellipse)
ax.text(2*np.pi/a-2.25/2, 0.36/2, "Ni", fontsize=20, ha='center', va='center', color='#7D2882')
ax.text(2*np.pi/a-0.25, 0.035, "Al", fontsize=14, ha='center', va='center', color='#F0AA23')

plt.savefig('plot/cylinder_60_Al_bg.png', dpi=600)

plt.show()
