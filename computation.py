import drawing
from classes import *
from typing import *
from generation import *
from statistics import mean
import math
import copy
import numpy


# Computes running mean of x over n values
def running_average(x, n: int = 3):
    ret = np.cumsum(x, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# Returns the energy level of the chain
# Epsilon is the energy associated with two H contacts
def calculate_energy(epsilon: float, lattice: ProteinLattice) -> float:
    # E = ε * f
    # E = Energy
    # ε = energy per H-H contact
    # f = Encounters
    f: int = 0

    for i in range(0, len(lattice.chain)):
        monomer = lattice.chain[i]

        # We only care about H monomers.
        if monomer.kind != MonomerKind.H:
            continue

        # Therefore we now are guaranteed to have an H monomer at chain[i].
        # Now we can check all of our neighbours using the lattice.
        for neighbour in lattice.get_neighbours(monomer.x, monomer.y):
            if neighbour.kind == MonomerKind.H:
                f += 1

    # Since we double count, we need to divide by 2.
    # Integer division is fine since we are guaranteed to be divisible by 2.
    f /= 2

    return -1.0 * epsilon * float(f)


# Returns the positions to check for a diff between previous and current points.
# Specifically for endpoint rotations.
def endpoints_rotate_lookup_table(diff_x: int, diff_y: int) -> List[Tuple[int, int]]:
    return {
        (0, +1): [(1, 0), (-1, 0)],  # If prev is above, return right/left
        (0, -1): [(1, 0), (-1, 0)],  # If prev is below, return right/left
        (+1, 0): [(0, 1), (0, -1)],  # If prev is right, return top/bottom
        (-1, 0): [(0, 1), (0, -1)],  # If prev is left, return top/bottom
    }[(diff_x, diff_y)]


# Returns position of potential kink jump, if any, based on neighbour positions.
# Returns List of coordinates with len = 1, or len = 0 if no potential kink jump.
def kink_jump_lookup_table(diff_prev_x: int,
                           diff_prev_y: int,
                           diff_next_x: int,
                           diff_next_y: int) -> List[Tuple[int, int]]:
    return {  # Only eight possible configurations can match, otherwise return empty list
        # Clock Wise
        ((0, 1), (1, 0)): [(1, 1)],  # prev above, next right
        ((1, 0), (0, -1)): [(1, -1)],  # prev right, next bottom
        ((0, -1), (-1, 0)): [(-1, -1)],  # prev bottom, next left
        ((-1, 0), (0, 1)): [(-1, 1)],  # prev left, next top

        # Counter Clock Wise
        ((0, 1), (-1, 0)): [(-1, 1)],  # prev above, next left
        ((-1, 0), (0, -1)): [(-1, -1)],  # prev left, next bottom
        ((0, -1), (1, 0)): [(1, -1)],  # prev bottom, next right
        ((1, 0), (0, 1)): [(1, 1)],  # prev right, next top
    }.get(((diff_prev_x, diff_prev_y), (diff_next_x, diff_next_y)), [])


# Tries to do a kink jump, returns whether it succeeded or not
# Might do an endpoint rotation if requested on endpoint in chain
# Modifies the given lattice!
def perform_kink_jump(idx_to_jump: int, lattice: ProteinLattice) -> bool:
    monomer = lattice.chain[idx_to_jump]
    if idx_to_jump == 0 or idx_to_jump == len(lattice.chain) - 1:
        # Endpoint rotate
        # Determine previous monomer
        if idx_to_jump == 0:
            prev_monomer = lattice.chain[1]
        else:
            prev_monomer = lattice.chain[-2]

        # Determine relative offsets from prev_monomer based on diff xy
        offsets = endpoints_rotate_lookup_table(prev_monomer.x - monomer.x, prev_monomer.y - monomer.y)
        for (x, y) in offsets:
            # Check if position is taken or not, and move if possible.
            if not lattice.has_monomer(prev_monomer.x + x, prev_monomer.y + y):
                lattice.move_monomer(idx_to_jump, prev_monomer.x + x, prev_monomer.y + y)
                return True
    else:
        # Determine if a kink jump can be performed
        prev_mon = lattice.chain[idx_to_jump - 1]
        next_mon = lattice.chain[idx_to_jump + 1]
        # Check the four configurations

        # Gather the position to check
        offsets = kink_jump_lookup_table(prev_mon.x - monomer.x,
                                         prev_mon.y - monomer.y,
                                         next_mon.x - monomer.x,
                                         next_mon.y - monomer.y)

        # Check if position is taken or not, perform kink jump if possible.
        for (x, y) in offsets:
            # List will always be either 0 or 1 in length.
            # (Using a list makes it easier to use than an explicit None check)
            if not lattice.has_monomer(monomer.x + x, monomer.y + y):
                lattice.move_monomer(idx_to_jump, monomer.x + x, monomer.y + y)
                return True

    # No kink jump possible. Returning False.
    return False


# Gathers all monomers that need to be rotated based on rotation point and part
def gather_rotated_monomers(lattice: ProteinLattice,
                            rotation_point_idx: int,
                            part: MonomerPart) -> List[Tuple[int, Monomer]]:
    monomers = []
    if part == MonomerPart.Left:
        # Get 0 - idx (exclusive idx)
        for i in range(0, rotation_point_idx):
            monomers.append((i, lattice.chain[i]))
    else:
        # Get idx - end (exclusive idx)
        for i in range(rotation_point_idx + 1, len(lattice.chain)):
            monomers.append((i, lattice.chain[i]))

    return monomers


# Tries to perform a pivot move given a rotation point, direction, which part to rotate and the input lattice.
# Returns true if rotation succeeded, false if not.
# Modifies the given lattice!
def perform_pivot(rotation_point_idx: int,
                  direction: Direction,
                  rotated_part: MonomerPart,
                  lattice: ProteinLattice) -> bool:
    rotation_monomer = lattice.chain[rotation_point_idx]

    rotated_part = gather_rotated_monomers(lattice, rotation_point_idx, rotated_part)

    # Exception for special case when rotating the endpoint and part to rotate is len == 0
    # In such a special case, nothing is done to the monomer, thus we would waste an iteration.
    # Therefore, we return False.
    if len(rotated_part) == 0:
        return False

    # Get the rotation angle
    angle = {
        Direction.ClockWise: 3 * math.pi / 2.0,
        Direction.CounterClockWise: math.pi / 2.0
    }[direction]

    new_positions = []
    for (mon_idx, monomer) in rotated_part:
        # Translate the monomer to 0,0
        shifted_x = monomer.x - rotation_monomer.x
        shifted_y = monomer.y - rotation_monomer.y

        # Rotate the monomer
        a1 = float(shifted_x) * math.cos(angle)
        a2 = float(shifted_y) * math.sin(angle)
        a3 = float(shifted_x) * math.sin(angle)
        a4 = float(shifted_y) * math.cos(angle)

        rotated_shifted_x = int(round(a1 - a2))
        rotated_shifted_y = int(round(a3 + a4))

        # Translate the monomer back to original position
        new_positions.append((mon_idx,
                              (rotated_shifted_x + rotation_monomer.x,
                               rotated_shifted_y + rotation_monomer.y)
                              ))

    # We need to check both the lattice if we can rotate.
    for (_, (x, y)) in new_positions:
        # Check for positions from the rotated part, these are always free after rotation.
        # (It effectively excludes the positions of the rotated part.)
        if any(elem.x == x and elem.y == y for (_, elem) in rotated_part):
            continue
        # Check the lattice
        if lattice.has_monomer(x, y):
            return False

    # Actually move the monomers to the new positions in the lattice/chain
    lattice.move_monomers(new_positions)

    return True


# Discards first {fraction} amount of elements from values.
def discard_fraction_of_array(values: List[float], fraction: float = 0.1) -> List[float]:
    slice_idx = int(math.ceil(len(values) * fraction))
    return values[slice_idx:len(values)]


# Computes the heat capacity using the input float array of energy levels
def compute_heat_capacity(values: List[float], temperature: float, boltzmann: float = 1.0,
                          discard: bool = True, discard_fraction: float = 0.1) -> float:
    if discard:
        values = discard_fraction_of_array(values, discard_fraction)

    mean_squared_energy = mean([e ** 2 for e in values])
    mean_energy = mean(values)

    return (mean_squared_energy - (mean_energy ** 2)) / (boltzmann * temperature)


# Performs the kink jump move as part of the main mmc loop.
# May fail, returns False in that case!
def mmc_attempt_kink_jump(lattice: ProteinLattice) -> bool:
    success = False
    # Perform kink jump

    # Keep track of attempted indices, since kink jumps and endpoint rotations are usually only possible on
    # a few positions in the chain. So therefore we keep track of this as a small optimization in order to avoid
    # re-attempting already attempted positions where we did a kink jump move.
    # If this list becomes empty, it means kink jump was not possible. The function will return False.
    # In such a situation the algorithm MUST perform a pivot first in order to get unstuck.
    possible_attempts = [e for e in range(0, len(lattice.chain))]
    while not success and len(possible_attempts) != 0:
        jump_idx = choices(possible_attempts)[0]
        success = perform_kink_jump(jump_idx, lattice)
        if not success:
            # Disallow attempted indices
            possible_attempts.remove(jump_idx)
    return success


# Performs the pivot move as part of the main mmc loop.
# In practice always succeeds so always should return True.
def mmc_perform_pivot(lattice: ProteinLattice) -> bool:
    success = False
    while not success:
        rotation_idx = choices(range(0, len(lattice.chain)))[0]
        direction = choice([0, 1])
        part = choice([0, 1])
        success = perform_pivot(rotation_idx, Direction(direction), MonomerPart(part), lattice)
    return success


# Returns the default protein
def mmc_initialize_default_protein(chain_length: int, hydrophobicity: float) -> ProteinLattice:
    # Generate the protein chain.
    seed(1234, 2)  # Set seed to fixed value for reproducibility of initial configuration.
    lattice = ProteinLattice(generate_protein(chain_length, hydrophobicity), hydrophobicity)

    return lattice


# Main function for performing the MMC simulation.
# Returns a tuple: ( (lowest_energy_lattice, lowest_energy), result_lattice, samples )
def mmc(temperature: float,  # Temperature parameter
        max_iterations: int,  # Iterations to do
        sampling_frequency: int,  # Frequency at which needs to be sampled
        lattice: ProteinLattice,  # Lattice that needs to be operated on
        epsilon: float = 1.0,  # Epsilon
        boltzmann: float = 1.0,  # Boltzmann weight
        draw_initial_conformation_plot: bool = False,  # Boolean indicating if initial conformation needs to be drawn
        draw_resulting_conformation_plot: bool = False,  # Boolean indicating if final conformation needs to be drawn
        # Keeps track of lowest lattice found. Causes noticable performance hit due to excess copying of memory.
        store_lowest_lattice: bool = False) -> Tuple[Tuple[ProteinLattice, float], ProteinLattice, MMCSamples]:

    # Draw the initial conformation or not
    if draw_initial_conformation_plot:
        drawing.draw_protein_conformation(lattice, temperature, lattice.hydrophobicity)

    energy_samples = []
    gyration_samples = []

    # Take initial samples
    energy = calculate_energy(epsilon, lattice)
    energy_samples.append(energy)
    gyration_samples.append(lattice.compute_gyration_radius())

    # Store which lattice is the lowest encountered so far.
    lowest_lattice = lattice
    lowest_lattice_energy: float = energy

    for iteration in range(0, max_iterations):
        # Choose operation
        operation_kind = choice([0, 1])
        if operation_kind == 0:
            # Perform kink jump / endpoint rotation.
            # In some rare cases this can fail, so we need to check for that.
            # In such situations there are no kink jump / endpoint rotations possible.
            # Therefore, opposed to the given sample pseudocode, I check this and perform a pivot instead.
            # This prevents the simulation from becoming stuck.
            success = mmc_attempt_kink_jump(lattice)
        else:
            # Perform pivot
            success = mmc_perform_pivot(lattice)

        # In certain rare cases a kink jump/endpoint_rotation is not possible,
        # so we need to perform a pivot instead.
        if not success and operation_kind == 0:
            mmc_perform_pivot(lattice)

        # We have successfully changed our chain here.
        new_energy = calculate_energy(epsilon, lattice)
        if new_energy < energy:
            energy = new_energy
            # Save the lattice as lowest using a deepcopy if this is new lowest energy lattice.
            if store_lowest_lattice and new_energy < lowest_lattice_energy:
                lowest_lattice = copy.deepcopy(lattice)
                lowest_lattice_energy = new_energy

        else:
            # boltzmann weight
            w: float = math.exp(- (float(new_energy) - float(energy)) / (float(boltzmann) * temperature))
            # Reject if w is smaller or equal to random value in 0-1
            if w > random():
                energy = new_energy
            else:
                lattice.undo_last_change()

        # Uncomment this line to print progress.
        # print('Iteration: {}/{} T: {}'.format(iteration + 1, max_iterations, temperature))

        # Sample the energy and gyration
        if (iteration + 1) % sampling_frequency == 0:
            energy_samples.append(energy)
            gyration_samples.append(lattice.compute_gyration_radius())

    # If enabled draw the resulting conformation plot
    if draw_resulting_conformation_plot:
        drawing.draw_protein_conformation(lattice, temperature, lattice.hydrophobicity)

    # Return values
    return (lowest_lattice, lowest_lattice_energy), lattice, MMCSamples(energy_samples, gyration_samples)
