from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from computation import *
from itertools import *

blue = np.array([65 / 256, 105 / 256, 225 / 256, 1])
orange = np.array([255 / 256, 165 / 256, 0 / 256, 1])
colormap = ListedColormap([blue, orange])


def plot_protein(lattice: ProteinLattice):
    #x_coords = [elem.x for elem in lattice.chain]
    #y_coords = [elem.y for elem in lattice.chain]
    #kinds = [int(elem.kind) for elem in lattice.chain]

    x_coords = [int(lattice.chain[elem.index].x)
                for elem in sorted([col for row in lattice.lattice for col in row], key=lambda e: e.index)
                if elem.index != -1]

    y_coords = [int(lattice.chain[elem.index].y)
                for elem in sorted([col for row in lattice.lattice for col in row], key=lambda e: e.index)
                if elem.index != -1]

    kinds = [int(lattice.chain[elem.index].kind)
             for elem in sorted([col for row in lattice.lattice for col in row], key=lambda e: e.index)
             if elem.index != -1]

    plt.title('Protein - N = {}, E = {}'.format(len(lattice.chain), calculate_energy(1.0, lattice)))
    plt.plot(x_coords, y_coords, 'k-', zorder=0)
    plt.scatter(x_coords, y_coords, c=kinds, cmap=colormap, zorder=3)
    plt.margins(0.1)
    plt.gca().set_aspect('equal')
    plt.show()
