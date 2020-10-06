from generation import *
from drawing import *


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    seed(1234, 2)  # Set seed to fixed value for reproducibility

    for i in range(0, 5):
        chain = generate_protein(25, 0.5)
        plot_protein(chain)
