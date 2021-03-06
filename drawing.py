from computation import *

blue = np.array([65 / 256, 105 / 256, 225 / 256, 1])
orange = np.array([255 / 256, 165 / 256, 0 / 256, 1])


# Compute next perfect square to determine grid size for histograms.
def next_perfect_square(N):
    next_n = math.floor(math.sqrt(N)) + 1
    return next_n * next_n


# Internal function for adjacency. Do not use.
def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


# Internal function used to set axis style of subplots.
def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)


# Draws violin plots vs. temperature.
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


# Plots histograms for energy/gyration vs. temperature.
# Plots them in a big figure.
def draw_histograms(temperatures: List[float],
                    values: List[List[float]],
                    xlabel: str,
                    ylabel: str,
                    title: str):
    cols = int(math.sqrt(next_perfect_square(len(temperatures))))
    rows = int(len(temperatures) / cols)
    print(len(temperatures))
    print(rows)
    print(cols)
    fig, ax = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(17, 23))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    for i in range(rows):
        for j in range(cols):
            if (i * rows) + j >= len(temperatures):
                break

            idx = (i * cols) + j
            axis = ax[i, j]

            axis.hist(values[idx], alpha=0.7, rwidth=0.85)
            axis.grid(axis='y', alpha=0.75)
            axis.set_title('T = {:.2f}'.format(temperatures[idx]))

    for axi in ax.flat:
        axi.set(xlabel=xlabel, ylabel=ylabel)

    for axi in ax.flat:
        axi.label_outer()

    fig.suptitle(title, y=0.99, size='large')
    fig.show()


# Plots the protein.
def draw_protein_conformation(lattice: ProteinLattice, temperature: float, hydrophobicity: float):
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
    # Sort results by temperature. min temp -> max temp
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

    # Draw distributions for energy vs temp
    if draw_energy_histograms_per_temp:
        draw_histograms(temperatures, [result_set[1] for result_set in results],
                        title='Energy distributions for different temperatures',
                        xlabel='Energy levels',
                        ylabel='Counts (relative)')

    # Draw distributions for gyration vs temp
    if draw_gyration_histograms_per_temp:
        draw_histograms(temperatures, [result_set[2] for result_set in results],
                        title='Gyration radius distributions for different temperatures',
                        xlabel='Gyration radii',
                        ylabel='Counts (relative)')

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
