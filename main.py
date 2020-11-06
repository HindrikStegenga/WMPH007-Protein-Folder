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
    hydrophobicity = 0.2
    lattice = mmc_initialize_default_protein(25, hydrophobicity)
    (lowest_lattice, lowest_energy, lowest_temp), lattice, results = perform_mmc_simulated_annealing(
        lattice,
        50,
        15000,
        2.0,
        randomize_seed=True,
        store_lowest_lattice=True)
    plot_protein(lowest_lattice, lowest_temp, hydrophobicity)
    print('Lowest protein: {}'.format(lowest_lattice.chain))
    draw_simulated_annealing_plots(lattice, results,
                                   draw_energy_histograms_per_temp=False,
                                   draw_gyration_histograms_per_temp=False)

    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
