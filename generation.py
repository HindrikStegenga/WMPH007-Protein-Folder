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


# Returns the count of a kind of monomer in the chain.
def get_kind_count(chain: List[Monomer], kind: MonomerKind) -> int:
    return sum(1 for c in chain if c.kind == kind)


# Returns a string containing H and P for each monomer
def get_chain_composition_string(chain: List[Monomer]) -> str:
    return_value = ''
    for m in chain:
        if m.kind == MonomerKind.H:
            return_value += 'H'
        else:
            return_value += 'P'
    return return_value


# Validates if a given position is occupied in the chain.
def has_monomer(x: int, y: int, chain: List[Monomer]) -> bool:
    for i in range(0, len(chain)):
        if chain[i].x == x and chain[i].y == y:
            return True
    return False


# Validates whether the given position is valid. (i.e. is not occupied)
# Also checks for dead ends in the lattice.
def is_valid_new_position(x: int, y: int, chain: List[Monomer], dead_chain: List[Monomer]) -> bool:
    if has_monomer(x, y, chain):
        return False
    if has_monomer(x, y, dead_chain):
        return False

    if is_dead_position(x, y, chain, dead_chain):
        return False

    return True


def is_dead_position(x: int, y: int, chain: List[Monomer], dead_chain: List[Monomer]) -> bool:
    # We now check for dead ends by walking through the chain searching for 4 neighbour points
    # If we can not move from this position to another, exclude it by returning False.
    if (
            has_monomer(x + 1, y, chain) and
            has_monomer(x - 1, y, chain) and
            has_monomer(x, y + 1, chain) and
            has_monomer(x, y - 1, chain)
    ):
        return True

    # We now check for dead ends by walking through the dead chain searching for 4 neighbour points
    # If we can not move from this position to another, exclude it by returning False.
    if (
            has_monomer(x + 1, y, dead_chain) and
            has_monomer(x - 1, y, dead_chain) and
            has_monomer(x, y + 1, dead_chain) and
            has_monomer(x, y - 1, dead_chain)
    ):
        return True
    return False


# Generates a random protein chain.
# length is the length of the chain.
# hydrophobicity is a fraction between 0 and 1 determining the relative amount of H monomers.
def generate_protein(length: int, hydrophobicity: float) -> List[Monomer]:
    # We keep track of all nodes we tried but ended up in a dead state.
    # This is so we can recursively track back until we find a valid path.
    dead_chain = []

    current_chain = [
        # Add initial monomer to make the algorithm simpler.
        Monomer(
            # Returns H or P depending on weight
            MonomerKind(choices([1, 2], [hydrophobicity, 1.0 - hydrophobicity])[0]),
            0,  # x coord
            0  # y coord
        )
    ]
    # We need to generate N - 1 additional monomers
    while True:
        # Generate direction with equal probability
        direction = choice([0, 1, 2, 3])
        (new_x, new_y) = generate_new_coords(direction, current_chain[-1].x, current_chain[-1].y)

        # Checks if we can add the monomer at the given position.
        if is_valid_new_position(new_x, new_y, current_chain, dead_chain):
            # Take a step if we can take it
            current_chain.append(
                Monomer(
                    # Returns H or P depending on weight
                    MonomerKind(choices([1, 2], [hydrophobicity, 1.0 - hydrophobicity])[0]),
                    new_x,
                    new_y
                )
            )
        else:
            dead_chain.append(Monomer(
                MonomerKind.H,
                new_x,
                new_y
            ))
            # Now we need to check if the previous position is dead or not, it might now be, in which case we need to
            # walk back
            prev = current_chain[-1]
            if is_dead_position(prev.x, prev.y, current_chain, dead_chain):
                current_chain.pop()
                dead_chain.append(prev)

        if length == len(current_chain):
            return current_chain
        continue
