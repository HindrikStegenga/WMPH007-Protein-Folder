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

# Calculates a new chain with a kink_jump_end_rot
# def calculate_new_chain_kink_jump_end_rot()
