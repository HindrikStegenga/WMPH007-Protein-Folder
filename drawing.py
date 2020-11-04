from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from computation import *

blue = np.array([65 / 256, 105 / 256, 225 / 256, 1])
orange = np.array([255 / 256, 165 / 256, 0 / 256, 1])


def plot_protein(lattice: ProteinLattice, temperature: float, hydrophobicity: float):
    plt.title('HP Protein, N = {}, E = {}, T = {}, H = {:.2}'.format(
        len(lattice.chain),
        calculate_energy(1.0, lattice),
        temperature,
        hydrophobicity
    ))

    plt.plot([elem.x for elem in lattice.chain],
             [elem.y for elem in lattice.chain],
             'k-', zorder=0)

    plt.scatter([elem.x for elem in lattice.chain if elem.kind == MonomerKind.H],
                [elem.y for elem in lattice.chain if elem.kind == MonomerKind.H],
                zorder=3, label='H', color=blue)

    plt.scatter([elem.x for elem in lattice.chain if elem.kind == MonomerKind.P],
                [elem.y for elem in lattice.chain if elem.kind == MonomerKind.P],
                zorder=3, label='P', color=orange)

    plt.margins(0.1)
    plt.ylabel('position y-coordinate')
    plt.xlabel('position x-coordinate')

    plt.legend()
    plt.gca().set_aspect('equal')
    plt.show()


# Draws the plot for energy vs. iterations
def draw_energy_iterations_plot(samples: List[float]):
    print(len(samples))
    plt.title('Energy vs. iterations')
    plt.plot(samples)
    plt.ylabel('Energy')
    plt.xlabel('Iterations (x100)')
    plt.show()