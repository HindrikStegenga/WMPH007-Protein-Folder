from classes import *
from typing import *
import math


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
        ((0, 1), (1, 0)): [(1, 1)],  # prev above, next right
        ((1, 0), (0, -1)): [(1, -1)],  # prev right, next bottom
        ((0, -1), (-1, 0)): [(-1, -1)],  # prev bottom, next left
        ((-1, 0), (0, 1)): [(-1, 1)]  # prev left, next top
    }.get(((diff_prev_x, diff_prev_y), (diff_next_x, diff_next_y)), [])


# Tries to do a kink jump, returns whether it succeeded or not
# Might do an endpoint rotation if requested on endpoint in chain
# Modifies the given lattice!
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
            if lattice.has_monomer(prev_monomer.x + x, prev_monomer.y + y):
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
            if lattice.has_monomer(monomer.x + x, monomer.y + x):
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
    # Get the rotation angle
    angle = {
        Direction.ClockWise: 3 * math.pi / 2.0,
        Direction.CounterClockWise: math.pi / 2.0
    }[direction]

    new_positions = []
    for (_, monomer) in rotated_part:
        # Shift all monomers
        shifted_x = monomer.x - rotation_monomer.x
        shifted_y = monomer.y - rotation_monomer.y

        rotated_shifted_x = int(round(float(shifted_x) * math.cos(angle)) - (float(shifted_y) * math.sin(angle)))
        rotated_shifted_y = int(round(float(shifted_x) * math.sin(angle)) + (float(shifted_y) * math.cos(angle)))

        # Shift monomers back
        new_positions.append((rotated_shifted_x + rotation_monomer.x,
                              rotated_shifted_y + rotation_monomer.y))

    # We need to check both the lattice if we can rotate.
    for (x, y) in new_positions:
        # Check for positions from the rotated part, these are always free after rotation.
        # (It effectively excludes the positions of the rotated part.)
        if any(elem.x == x and elem.y == y for (_, elem) in rotated_part):
            continue
        # Check the lattice
        if lattice.has_monomer(x, y):
            return False

    # Actually move the monomers to the new positions in the lattice/chain
    for idx in range(0, len(rotated_part)):
        (idx_in_chain, _) = rotated_part[idx]
        (x, y) = new_positions[idx]
        lattice.move_monomer(idx_in_chain, x, y)

    return True
