from generation import *
from drawing import *


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    seed(1234, 2)  # Set seed to fixed value for reproducibility

    for i in range(0, 1):
        chain = generate_protein(25, 0.5)
        lattice = ProteinLattice(chain)
        print(perform_kink_jump(9, lattice))
        plot_protein(lattice)
        # for i in chain:
            # print(lattice.get_by_coordinate(i.x, i.y))
