from classes import *
from typing import *


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
    return {  # Only four possible configurations can match, otherwise return empty list
        ((0,  1), (1,  0)): [(1, 1)],    # prev above, next right
        ((1,  0), (0, -1)): [(1, -1)],   # prev right, next bottom
        ((0, -1), (-1, 0)): [(-1, -1)],  # prev bottom, next left
        ((-1, 0), (0,  1)): [(-1, 1)]    # prev left, next top
    }.get(((diff_prev_x, diff_prev_y), (diff_next_x, diff_next_y)), [])


# Tries to do a kink jump, returns whether it succeeded or not
# Might do an endpoint rotation if requested on endpoint in chain
# Modifies the inserted lattice!
# TODO: Implement random selection of endpoint rotation.
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
            if lattice.lattice[prev_monomer.x + x][prev_monomer.y + y].index == -1:
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
            if lattice.lattice[monomer.x + x][monomer.y + x].index == -1:
                lattice.move_monomer(idx_to_jump, monomer.x + x, monomer.y + y)
                return True

    # No kink jump possible. Returning False.
    return False


# Tries to perform a pivot move
def perform_pivot(endpoint: Monomer, rotation: Rotation, lattice: ProteinLattice) -> bool:
    return True

# Calculates a new chain with a kink_jump
# def calculate_new_chain_kink_jump_end_rot()
