from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from computation import *

blue = np.array([65 / 256, 105 / 256, 225 / 256, 1])
orange = np.array([255 / 256, 165 / 256, 0 / 256, 1])


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)


def draw_violin_plot_over_temp(title: str,
                               ylabel: str,
                               values: List[List[float]],
                               temperatures: List[float]):
    fig, axis = plt.subplots(nrows=1, ncols=1)
    axis.set_title(title)
    parts = axis.violinplot(
        values, showmeans=False, showmedians=False,
        showextrema=False)

    for pc in parts['bodies']:
        pc.set_facecolor('#D43F3A')
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    quartile1, medians, quartile3 = np.percentile(values, [25, 50, 75], axis=1)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(values, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    axis.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    axis.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    axis.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
    set_axis_style(axis, ['{:.1f}'.format(e) for e in temperatures])
    axis.set_ylabel(ylabel)
    axis.set_xlabel('Temperature (ε/kB)')
    plt.show()




def plot_protein(lattice: ProteinLattice, temperature: float, hydrophobicity: float):
    plt.title('HP Protein, N = {}, E = {:.2f}, T = {:.2f}, H = {:.2f}'.format(
        len(lattice.chain),
        calculate_energy(1.0, lattice),
        temperature,
        hydrophobicity
    ))

    plt.plot([elem.x for elem in lattice.chain],
             [elem.y for elem in lattice.chain],
             'k-', zorder=0)

    plt.scatter([elem.x for elem in lattice.chain if elem.kind == MonomerKind.H],
                [elem.y for elem in lattice.chain if elem.kind == MonomerKind.H],
                zorder=3, label='H', color=blue)

    plt.scatter([elem.x for elem in lattice.chain if elem.kind == MonomerKind.P],
                [elem.y for elem in lattice.chain if elem.kind == MonomerKind.P],
                zorder=3, label='P', color=orange)

    plt.margins(0.1)
    plt.ylabel('position y-coordinate')
    plt.xlabel('position x-coordinate')

    plt.legend()
    plt.gca().set_aspect('equal')
    plt.show()


# Draws the plot for energy vs. iterations
def draw_energy_iterations_plot(samples: List[float]):
    print(len(samples))
    plt.title('Energy vs. iterations')
    plt.plot(samples)
    plt.ylabel('Energy')
    plt.xlabel('Iterations (x100)')
    plt.show()


# Draws plots for simulated annealing
def draw_simulated_annealing_plots(lattice: ProteinLattice,
                                   results: List[Tuple[float, List[float], List[float]]],
                                   draw_energy_histograms_per_temp: bool = False,
                                   draw_gyration_histograms_per_temp: bool = False):
    results.sort(key=lambda x: x[0])
    temperatures = [elem[0] for elem in results]

    # Compute averages
    avg_energy_per_temp: List[float] = []
    avg_gyration_per_temp: List[float] = []
    for result_set in results:
        avg_energy = mean(result_set[1])
        avg_gyration = mean(result_set[2])

        avg_energy_per_temp.append(avg_energy)
        avg_gyration_per_temp.append(avg_gyration)

    # Draw avg energy vs temp
    plt.plot(temperatures, avg_energy_per_temp)
    plt.title('Average Energy vs. Temperature')
    plt.xlabel('Temperature (ε/kB)')
    plt.ylabel('Average energy')
    plt.show()

    # Draw avg gyration vs temp
    plt.plot(temperatures, avg_gyration_per_temp)
    plt.title('Average gyration vs. Temperature')
    plt.xlabel('Temperature (ε/kB)')
    plt.ylabel('Average gyration')
    plt.show()

    # Draw distributions for gyration vs temp
    if draw_energy_histograms_per_temp:
        for result_set in results:
            plt.hist(result_set[1], alpha=0.7, rwidth=0.85)
            plt.grid(axis='y', alpha=0.75)
            plt.xlabel('Energy radius')
            plt.ylabel('Counts')
            plt.title('Energy distribution at T = {:.2f}'.format(result_set[0]))
            plt.show()

    # Draw distributions for gyration vs temp
    if draw_gyration_histograms_per_temp:
        for result_set in results:
            plt.hist(result_set[2], alpha=0.7, rwidth=0.85)
            plt.grid(axis='y', alpha=0.75)
            plt.xlabel('Gyration radius')
            plt.ylabel('Counts')
            plt.title('Gyration radius distribution at T = {:.2f}'.format(result_set[0]))
            plt.show()

    # Draw violinplot for energy distributions
    draw_violin_plot_over_temp('Energy distributions per temperature',
                               'Energy level',
                               [elem[1] for elem in results],
                               temperatures)

    # Draw violinplot for gyration distributions
    draw_violin_plot_over_temp('Gyration distributions per temperature',
                               'Gyration radius',
                               [elem[2] for elem in results],
                               temperatures)

    # Compute heat capacity and draw plot vs temperature
    heat_capacity_per_temp: List[float] = [compute_heat_capacity(elem[1], elem[0])
                                           for elem in results]
    plt.plot(temperatures, heat_capacity_per_temp)
    plt.xlabel('Temperature (ε/kB)')
    plt.ylabel('Heat capacity')
    plt.title('Heat capacity vs. Temperature')
    plt.show()
