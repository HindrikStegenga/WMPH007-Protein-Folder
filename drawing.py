from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from computation import *

blue = np.array([65 / 256, 105 / 256, 225 / 256, 1])
orange = np.array([255 / 256, 165 / 256, 0 / 256, 1])
colormap = ListedColormap([blue, orange])


def plot_protein(lattice: ProteinLattice):
    x_coords = [elem.x for elem in lattice.chain]
    y_coords = [elem.y for elem in lattice.chain]
    kinds = [int(elem.kind) for elem in lattice.chain]

    plt.title('Protein - N = {}, E = {}'.format(len(lattice.chain), calculate_energy(1.0, lattice)))
    plt.plot(x_coords, y_coords, 'k-', zorder=0)
    plt.scatter(x_coords, y_coords, c=kinds, cmap=colormap, zorder=3)
    plt.margins(0.1)
    plt.gca().set_aspect('equal')
    plt.show()

