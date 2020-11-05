from debug_functions import check_random_walk
from benchmarking import *
import multiprocessing
import statistics
from collections import Counter


# Function that performs the MMC simulations, is used only as a way to execute simulations in parallel.
# This is done by spawning a separate process using the multiprocessing library
# This massively speeds up our total time required for simulations.
def perform_mmc(temp: float, out_list: List[AMMCSamples], sample_count: int = 1):
    c_samples = []
    samples = []
    for i in range(0, sample_count):
        _, samples = mmc(temp, 25, 20000, 100, 0.5, 1.0, 1.0, True)
        c_samples.append(compute_heat_capacity(samples.energy, temp))
        print('Computed MMC for T: {:.4} - {}/{}'.format(temp, i + 1, sample_count))

    out_list.append(AMMCSamples(samples.energy, samples.gyration_radius, mean(c_samples), temp))


# This function is ran for the simulated annealing procedure
def perform_mmc_simulated_annealing():
    # Total steps sampled from the temperature range
    mmc_steps = 50
    # Samples to take per mmc temperature step
    samples_per_step = 1
    # Maximum temperature to use.
    max_temp = 2.0
    # Minimum temperature to use.
    min_temp = 0.0

    # Temperature step per mmc step
    temp_step = (max_temp - min_temp) / mmc_steps
    # Job list for multi processing
    jobs = []

    # Set up a return-value array using shared memory between multiple subprocesses
    # This is where results of simulations are stored.
    # The individual MMC simulations are independent, so we can easily parallelize it.
    output_list: [AMMCSamples] = multiprocessing.Manager().list()

    # Set up process pool
    pool = multiprocessing.Pool(multiprocessing.cpu_count())

    # Iterate of the temperature range from min temp to max temp
    for i in range(0, mmc_steps):
        # Compute temperature, since temperatures are separately simulated it doesn't matter which order we start in.
        temperature = min_temp + float(i * temp_step)
        # Spawn job as a process to be submitted onto the process pool
        job = pool.apply_async(perform_mmc, (temperature, output_list, samples_per_step))
        jobs.append(job)

    # Execute the jobs
    pool.close()
    # Wait for the jobs
    pool.join()
    # Sort the results since order of threads is undefined
    output_list = sorted(output_list, key=lambda x: x.temperature)
    # print([e.temperature for e in output_list])
    # Gather the samples for each value
    total_e_samples: [float] = []
    for samples in output_list:
        discarded_g_radii = discard_fraction_of_array(samples.energy)
        total_e_samples += discarded_g_radii

    total_g_samples: [float] = []
    for samples in output_list:
        discarded_g_radii = discard_fraction_of_array(samples.gyration_radius)
        total_g_samples += discarded_g_radii

    c_samples: [float] = [e.heat_capacity for e in output_list]
    t_samples: [float] = [e.temperature for e in output_list]

    # Plot the results
    plt.title('Heat capacity vs. temperature')
    plt.plot(t_samples, c_samples)
    plt.ylabel('Heat capacity')
    plt.xlabel('Temperature (ε/kB)')
    plt.show()
    return
    energy_levels = set()
    temperatures = [elem.temperature for elem in output_list]
    for elem in output_list:
        for energy_level in elem.energy:
            energy_levels.add(energy_level)

    # plt.hist(total_e_samples, alpha=0.7, rwidth=0.85)

    for smmc in output_list:
        #     plt.hist(smmc.energy, bins=len(energy_levels), alpha=0.7, rwidth=0.85)
        #     plt.grid(axis='y', alpha=0.75)
        #     plt.xlabel('Energy')
        #     plt.ylabel('Counts')
        #     plt.title('Energy distribution at T = {}'.format(smmc.temperature))
        #     plt.show()

        plt.hist(smmc.gyration_radius, alpha=0.7, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Gyration radius')
        plt.ylabel('Counts')
        plt.title('Gyration radius distribution at T = {}'.format(smmc.temperature))
        plt.show()

    # plt.violinplot([elem.energy for elem in output_list], showmeans=True)
    # plt.hist(total_e_samples, alpha=0.7, rwidth=0.85)
    # plt.grid(axis='y', alpha=0.75)
    # plt.xlabel('Temperature')
    # plt.ylabel('Counts')
    # plt.title('Energy distribution')
    # plt.xticks(np.arange(min_temp, max_temp, step=0.1))
    # plt.show()

    # plt.hist(total_g_samples, alpha=0.7, rwidth=0.85)
    # plt.grid(axis='y', alpha=0.75)
    # plt.xlabel('Gyration radius')
    # plt.ylabel('Counts')
    # plt.title('Gyration radius distribution')
    # plt.show()

    # plt.hist2d([elem.temperature for elem in output_list], [elem.energy for elem in output_list])
    # plt.show()


def perform_mmc_simulated_annealing2(
        lattice: ProteinLattice,
        temperature_steps: int,  # Total temperature steps
        mmc_iterations_per_step: int,  # MMC iterations_per_temperature.
        max_temp: float,  # Maximum temperature to use.
        min_temp: float = 0.0,  # Minimum temperature to use.
        sampling_frequency: int = 100,
        epsilon: float = 1.0,
        boltzmann: float = 1.0,
        randomize_seed: bool = True, ) -> Tuple[ProteinLattice, List[Tuple[float, List[float], List[float]]]]:
    # (final_temp, energy[], gyration[])
    results: List[Tuple[float, List[float], List[float]]] = []

    # Set seed for MMC
    if randomize_seed:
        seed()

    # Temperature step per mmc step
    for iteration in range(0, temperature_steps):
        # Compute temperature
        temperature = max_temp - (((max_temp - min_temp) / temperature_steps) * iteration)
        print('Annealing at T: {:.2f}, {}/{}...'.format(temperature, iteration + 1, temperature_steps))

        # Perform mmc at the given temperature
        _, samples = mmc(temperature, mmc_iterations_per_step, sampling_frequency, lattice,
                         draw_initial_conformation_plot=(iteration == 0),
                         draw_resulting_conformation_plot=(iteration == temperature_steps - 1),
                         epsilon=epsilon,
                         boltzmann=boltzmann)

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

    return lattice, results
