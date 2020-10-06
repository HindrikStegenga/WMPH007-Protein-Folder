from typing import *
from classes import *
from random import *


# Directions as per random number generator:
# 0 => up       (x, y+1)
# 1 => right    (x+1, y)
# 2 => left     (x-1, y)
# 3 => down     (x, y-1)

# ProteinKind as per random number generator:
# 1 => H
# 2 => P

# Generates new coords based on a direction integer mapping and the old coords.
# direction MUST be part of {0,1,2,3}
def generate_new_coords(direction: int, old_x: int, old_y: int) -> (int, int):
    return {
        0: (old_x, old_y + 1),
        1: (old_x + 1, old_y),
        2: (old_x - 1, old_y),
        3: (old_x, old_y - 1)
    }[direction]


# Validates if a given position is occupied in the chain.
def has_monomer(x: int, y: int, chain: List[ProteinMonomer]) -> bool:
    for i in range(0, len(chain)):
        if chain[i].x == x and chain[i].y == y:
            return True
    return False


# Validates whether the given position is valid. (i.e. is not occupied)
# Also checks for dead ends in the lattice.
def is_valid_new_position(x: int, y: int, chain: List[ProteinMonomer]) -> bool:
    if has_monomer(x, y, chain):
        return False

    # We now check for dead ends by walking through the chain searching for 4 neighbour points
    # If we can not move from this position to another, exclude it by returning False.
    if (
            has_monomer(x + 1, y, chain) and
            has_monomer(x - 1, y, chain) and
            has_monomer(x, y + 1, chain) and
            has_monomer(x, y - 1, chain)
    ):
        return False

    return True


# Generates a random protein chain.
# length is the length of the chain.
# hydrophobicity is a fraction between 0 and 1 determining the relative amount of H monomers.
def generate_protein(length: int, hydrophobicity: float) -> List[ProteinMonomer]:
    current_chain = [
        # Add initial monomer to make the algorithm simpler.
        ProteinMonomer(
            # Returns H or P depending on weight
            ProteinKind(choices([1, 2], [hydrophobicity, 1.0 - hydrophobicity])[0]),
            0,  # x coord
            0  # y coord
        )
    ]
    # We need to generate N - 1 additional monomers
    for i in range(0, length - 1):
        while True:
            # Generate direction with equal probability
            direction = choice([0, 1, 2, 3])
            (new_x, new_y) = generate_new_coords(direction, current_chain[-1].x, current_chain[-1].y)
            if is_valid_new_position(new_x, new_y, current_chain):
                current_chain.append(
                    ProteinMonomer(
                        # Returns H or P depending on weight
                        ProteinKind(choices([1, 2], [hydrophobicity, 1.0 - hydrophobicity])[0]),
                        new_x,
                        new_y
                    )
                )
                break
            else:
                continue

    return current_chain
