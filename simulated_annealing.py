from debug_functions import check_random_walk
from benchmarking import *
import multiprocessing
import statistics
from collections import Counter


def perform_mmc_simulated_annealing(
        lattice: ProteinLattice,
        temperature_steps: int,  # Total temperature steps
        mmc_iterations_per_step: int,  # MMC iterations_per_temperature.
        max_temp: float,  # Maximum temperature to use.
        min_temp: float = 0.0,  # Minimum temperature to use.
        sampling_frequency: int = 100,
        epsilon: float = 1.0,
        boltzmann: float = 1.0,
        randomize_seed: bool = True,
        store_lowest_lattice: bool = False) -> Tuple[Tuple[ProteinLattice, float, float],
                                                     ProteinLattice, List[Tuple[float,
                                                                                List[float],
                                                                                List[float]]]]:
    # (final_temp, energy[], gyration[])
    results: List[Tuple[float, List[float], List[float]]] = []

    # Set seed for MMC
    if randomize_seed:
        seed()

    lowest_lattice = lattice
    lowest_lattice_energy: float = calculate_energy(epsilon, lattice)
    lowest_temp: float = max_temp

    # Temperature step per mmc step
    for iteration in range(0, temperature_steps):
        # Compute temperature
        temperature = max_temp - (((max_temp - min_temp) / temperature_steps) * iteration)
        print('Annealing at T: {:.2f}, {}/{}...'.format(temperature, iteration + 1, temperature_steps))

        # Perform mmc at the given temperature
        (lowest, lowest_energy), _, samples = mmc(temperature, mmc_iterations_per_step, sampling_frequency, lattice,
                                                  draw_initial_conformation_plot=(iteration == 0),
                                                  draw_resulting_conformation_plot=(iteration == temperature_steps - 1),
                                                  epsilon=epsilon,
                                                  boltzmann=boltzmann,
                                                  store_lowest_lattice=store_lowest_lattice)

        if store_lowest_lattice and lowest_energy < lowest_lattice_energy:
            lowest_lattice = copy.deepcopy(lowest)
            lowest_lattice_energy = lowest_energy
            lowest_temp = temperature

        # Generate results, discarding first 10%
        # (final_temp, energy[], gyration[])
        results.append((temperature, discard_fraction_of_array(samples.energy),
                        discard_fraction_of_array(samples.gyration_radius)))

    print('Annealing at T: {:.2f}, {}/{}... done.'.format(min_temp, temperature_steps, temperature_steps))
    print('Final energy: {}'.format(calculate_energy(epsilon, lattice)))
    all_energy = [val
                  for result_set in results
                  for val in result_set[1]]

    all_gyration = [val
                    for result_set in results
                    for val in result_set[2]]

    print('Lowest energy state found: {:.2f}'.format(min(all_energy)))
    print('Mean energy state: {:.2f}'.format(mean(all_energy)))

    print('Lowest gyration radius found: {:.2f}'.format(min(all_gyration)))
    print('Mean gyration radius: {:.2f}'.format(mean(all_gyration)))

    print('Resulting protein: ')
    print(lattice.chain)

    return (lowest_lattice, lowest_lattice_energy, lowest_temp), lattice, results
