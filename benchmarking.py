from drawing import *


# This function is run for the benchmarking procedure
def perform_mmc_benchmarking():
    iterations = 50000
    length = 25
    temperature = 0.25
    sample_frequency = 100
    averaging = 3

    lattice02 = mmc_initialize_default_protein(length, 0.2)
    lattice05 = mmc_initialize_default_protein(length, 0.5)
    lattice08 = mmc_initialize_default_protein(length, 0.8)

    seed()  # Set fresh seed
    _, samples = mmc(temperature,  # Temperature
                     iterations,  # Total amount of iterations
                     sample_frequency,  # Iterations at which to sample
                     lattice02,
                     draw_initial_conformation_plot=False,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(running_average(samples.energy, averaging))

    seed()  # Set fresh seed
    _, samples = mmc(temperature,  # Temperature
                     iterations,  # Total amount of iterations
                     sample_frequency,  # Iterations at which to sample
                     lattice05,
                     draw_initial_conformation_plot=False,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(running_average(samples.energy, averaging))

    seed()  # Set fresh seed
    _, samples = mmc(temperature,  # Temperature
                     iterations,  # Total amount of iterations
                     sample_frequency,  # Iterations at which to sample
                     lattice08,
                     draw_initial_conformation_plot=False,
                     draw_resulting_conformation_plot=True)
    draw_energy_iterations_plot(running_average(samples.energy, averaging))
