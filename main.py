from debug_functions import check_random_walk
from drawing import *
import multiprocessing


# Function that performs the MMC simulations, is used only as a way to execute simulations in parallel.
# This is done by spawning a separate process using the multiprocessing library
# This massively speeds up our total time required for simulations.
def perform_mmc(temp: float, out_list: List[Tuple[float, float]]):
    _, es = mmc(temp, 25, 40000, 100, 0.5, 1.0, 1.0, True)
    out_list.append((temp, compute_heat_capacity(es, temp)))
    print('Computed MMC for T: {}'.format(temp))


# Draws the plot for energy vs. iterations
def draw_energy_iterations_plot(samples: List[float]):
    plt.title('Energy vs. iterations')
    plt.plot(samples)
    plt.ylabel('Energy (ε/kB)')
    plt.xlabel('Iterations (x100)')
    plt.show()


# This function is run for the benchmarking procedure
def perform_mmc_benchmarking():
    _, samples = mmc(0.2,  # Temperature
                     25,  # Length of the chain
                     35000,  # Total amount of iterations
                     100,  # Iterations at which to sample
                     0.2,  # Hydrophobic fraction
                     draw_initial_conformation_plot=True,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(samples)
    _, samples = mmc(0.2,  # Temperature
                     25,  # Length of the chain
                     35000,  # Total amount of iterations
                     100,  # Iterations at which to sample
                     0.5,  # Hydrophobic fraction
                     draw_initial_conformation_plot=True,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(samples)
    _, samples = mmc(0.2,  # Temperature
                     25,  # Length of the chain
                     35000,  # Total amount of iterations
                     100,  # Iterations at which to sample
                     0.8,  # Hydrophobic fraction
                     draw_initial_conformation_plot=True,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(samples)


# This function is ran for the simulated annealing procedure
def perform_mmc_simulated_annealing():
    # Total steps sampled from the temperature range
    mmc_steps = 200
    # Maximum temperature to set.
    max_temp = 3.0
    jobs = []

    # Set up a return-value array using shared memory between multiple subprocesses
    # This is where results of simulations are stored.
    # The individual MMC simulations are independent, so we can easily parallelize it.
    output_list = multiprocessing.Manager().list()

    # Set up process pool
    pool = multiprocessing.Pool(multiprocessing.cpu_count())

    # Iterate of the temperature range downwards from max_temp to 0
    for i in range(mmc_steps, 0, -1):
        # Compute temperature
        temperature = float(i * max_temp) / float(mmc_steps)
        # Spawn job as a process to be submitted onto the process pool
        job = pool.apply_async(perform_mmc, (temperature, output_list))
        jobs.append(job)

    # Execute the jobs
    pool.close()
    # Wait for the jobs
    pool.join()
    # Sort the results since order of threads is undefined
    output_list = sorted(output_list, key=lambda x: x[0])

    # Gather the samples for each value
    c_samples = [c for t, c in output_list]
    t_samples = [t for t, c in output_list]

    # Plot the results
    plt.title('Heat capacity vs. temperature')
    plt.plot(t_samples, c_samples)
    plt.ylabel('heat capacity')
    plt.xlabel('temperature')
    plt.show()


# Main function of the program
def main():
    # Enable following two lines to run the benchmarking procedure for 3 fixed temperatures.
    # perform_mmc_benchmarking()
    # return

    # Enable following two lines to run the simulated annealing procedure.
    perform_mmc_simulated_annealing()
    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
