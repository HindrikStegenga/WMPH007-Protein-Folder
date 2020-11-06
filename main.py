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
    hydrophobicity = 0.5
    lattice = mmc_initialize_default_protein(25, hydrophobicity)
    (lowest_lattice, lowest_energy, lowest_temp), lattice, results = perform_mmc_simulated_annealing(
        lattice,
        25,  # Temperature steps
        15000,  # MMC Iterations per temperature step
        2.0,  # Max temperature
        min_temp= 0.0,  # Min temperature
        randomize_seed=True,
        store_lowest_lattice=False)  # Set this only to true if necessary, slows down the algorithm significantly
    # Draw lowest protein
    draw_protein_conformation(lowest_lattice, lowest_temp, hydrophobicity)
    print('Lowest protein: {}'.format(lowest_lattice.chain))
    # Draw plots for results of annealing
    draw_simulated_annealing_plots(lattice, results,
                                   draw_energy_histograms_per_temp=False,
                                   draw_gyration_histograms_per_temp=False)

    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
