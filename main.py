from generation import *
from drawing import *


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    seed(1234, 2)  # Set seed to fixed value for reproducibility
    mmc(1, 25, 1000, 100, 0.5)
    #chain = generate_protein(25, 0.5)
    #lattice = ProteinLattice(chain)
    #print(perform_kink_jump(0, lattice))
    #plot_protein(lattice)
