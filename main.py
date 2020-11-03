from debug_functions import check_random_walk
from drawing import *
import multiprocessing


# Function that performs the MMC simulations, is used only as a way to execute simulations in parallel.
# This is done by spawning a separate process using the multiprocessing library
# This massively speeds up our total time required for simulations.
def perform_mmc(temp: float, out_list: List[Tuple[float, float]]):
    _, es = mmc(temp, 25, 30000, 100, 0.5, 1.0, 1.0, True)
    out_list.append((temp, compute_heat_capacity(es, temp)))
    print('MMC: {}'.format(temp))


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


# Main function of the program
def main():
    perform_mmc_benchmarking()
    return

    mmc_steps = 50
    jobs = []
    manager = multiprocessing.Manager()
    output_list = manager.list()
    for i in range(mmc_steps, 0, -1):
        temperature = float(i * 2) / float(mmc_steps)
        process = multiprocessing.Process(target=perform_mmc, args=(temperature, output_list))
        jobs.append(process)

    for job in jobs:
        job.start()
    for job in jobs:
        job.join()

    output_list = sorted(output_list, key=lambda x: x[0])
    print(output_list)
    c_samples = [c for t, c in output_list]
    t_samples = [t for t, c in output_list]

    plt.title('Energy vs. iterations')
    plt.plot(t_samples, c_samples)
    plt.ylabel('heat capacity')
    plt.xlabel('temperature')
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
