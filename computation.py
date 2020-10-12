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


# Tries to do a kink jump, returns whether it succeeded or not
# Might do an endpoint rotation if requested on endpoint in chain
# Modifies the inserted lattice!
def perform_kink_jump(idx_to_jump: int, lattice: ProteinLattice) -> bool:
    monomer = lattice.chain[idx_to_jump]
    if idx_to_jump == 0 or idx_to_jump == len(lattice.chain) - 1:
        # Endpoint rotate

        # Determine previous monomer
        if idx_to_jump == 0:
            prev_monomer = lattice.chain[1]
        else:
            prev_monomer = lattice.chain[-2]

        # Determine if endpoint rotate can be performed
        if ((prev_monomer.x == monomer.x and prev_monomer.y == monomer.y + 1) or
                (prev_monomer.x == monomer.x and prev_monomer.y == monomer.y - 1)):
            # Prev = top/bottom, try rotate right, then left. (clockwise)
            if lattice.lattice[prev_monomer.x + 1][prev_monomer.y].index == -1:
                # Right
                lattice.move_monomer(idx_to_jump, prev_monomer.x + 1, prev_monomer.y)
                return True
            if lattice.lattice[prev_monomer.x - 1][prev_monomer.y].index == -1:
                # Left
                lattice.move_monomer(idx_to_jump, prev_monomer.x - 1, prev_monomer.y)
                return True

        if ((prev_monomer.x + 1 == monomer.x and prev_monomer.y == monomer.y) or
                (prev_monomer.x - 1 == monomer.x and prev_monomer.y == monomer.y)):
            # Prev = left/right, try rotate top, then bottom. (clockwise)
            if lattice.lattice[prev_monomer.x][prev_monomer.y + 1].index == -1:
                # Top
                lattice.move_monomer(idx_to_jump, prev_monomer.x, prev_monomer.y + 1)
                return True
            if lattice.lattice[prev_monomer.x][prev_monomer.y - 1].index == -1:
                # Bottom
                lattice.move_monomer(idx_to_jump, prev_monomer.x - 1, prev_monomer.y)
                return True

    else:
        # Determine if a kink jump can be performed
        prev_mon = lattice.chain[idx_to_jump - 1]
        next_mon = lattice.chain[idx_to_jump + 1]
        # Check the four configurations
        # Check
        if (prev_mon.x == monomer.x and
                prev_mon.y == monomer.y + 1 and
                next_mon.x == monomer.x + 1 and
                next_mon.y == monomer.y):
            # prev = top, next is right
            # Check x + 1, y + 1
            if lattice.lattice[monomer.x + 1][monomer.y + 1].index == -1:
                # Hooray we can do kink jump here!
                lattice.move_monomer(idx_to_jump, monomer.x + 1, monomer.y + 1)
                return True

        if (prev_mon.x == monomer.x + 1 and
                prev_mon.y == monomer.y and
                next_mon.x == monomer.x and
                next_mon.y == monomer.y - 1):
            # prev = right, next is bottom
            # Check x + 1, y - 1
            if lattice.lattice[monomer.x + 1][monomer.y - 1].index == -1:
                # Hooray we can do kink jump here!
                lattice.move_monomer(idx_to_jump, monomer.x + 1, monomer.y - 1)
                return True

        if (prev_mon.x == monomer.x and
                prev_mon.y == monomer.y - 1 and
                next_mon.x == monomer.x - 1 and
                next_mon.y == monomer.y):
            # prev = bottom, next is left
            # Check x - 1, y - 1
            if lattice.lattice[monomer.x - 1][monomer.y - 1].index == -1:
                # Hooray we can do kink jump here!
                lattice.move_monomer(idx_to_jump, monomer.x - 1, monomer.y - 1)
                return True

        if (prev_mon.x == monomer.x - 1 and
                prev_mon.y == monomer.y and
                next_mon.x == monomer.x and
                next_mon.y == monomer.y + 1):
            # prev = left, next is top
            # Check x - 1, y + 1
            if lattice.lattice[monomer.x - 1][monomer.y + 1].index == -1:
                # Hooray we can do kink jump here!
                lattice.move_monomer(idx_to_jump, monomer.x - 1, monomer.y + 1)
                return True

    return False


# Tries to perform a pivot move
def perform_pivot(endpoint: Monomer, rotation: Rotation, lattice: ProteinLattice) -> bool:
    return True

# Calculates a new chain with a kink_jump
# def calculate_new_chain_kink_jump_end_rot()
