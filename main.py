from debug_functions import check_random_walk
from benchmarking import *
from simulated_annealing import *


# Main function of the program
def main():
    # Enable following two lines to run the random walk checking procedure
    # check_random_walk()
    # return

    # Enable following two lines to run the benchmarking procedure for 3 fixed temperatures.
    # perform_mmc_benchmarking()
    # return

    # Enable following lines to run the simulated annealing procedure.
    lattice = mmc_initialize_default_protein(25, 0.5)
    lattice, results = perform_mmc_simulated_annealing2(lattice, 20, 15000, 2.0)
    draw_simulated_annealing_plots(lattice, results)

    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
