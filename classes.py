from matplotlib import pyplot as plt
import numpy as np
from enum import IntEnum
from typing import *
import math


class Direction(IntEnum):
    ClockWise = 0
    CounterClockWise = 1


class MonomerPart(IntEnum):
    Left = 0
    Right = 1


class MonomerRecordValue(IntEnum):
    NONE = 0
    H = 1
    P = 2


class MonomerMoveRecord:
    def __init__(self, index: int, old: Tuple[int, int], new: Tuple[int, int]):
        self.index = index
        self.old = old
        self.new = new


class MonomerRecord:
    def __init__(self, value: MonomerRecordValue, index_in_chain: int):
        self.value: MonomerRecordValue = value
        self.index: int = index_in_chain


class MonomerKind(IntEnum):
    H = 1
    P = 2

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return {
            1: 'H',
            2: 'P'
        }[self]


class Monomer:
    def __init__(self, kind: MonomerKind, x: int, y: int):
        self.kind = kind
        self.x: int = x
        self.y: int = y

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '(({},{}) {})'.format(self.x, self.y, self.kind.__str__())


def calculate_chain_dimensions(chain: List[Monomer]) -> (int, int):
    min_x = min([elem.x for elem in chain])
    min_y = min([elem.y for elem in chain])
    max_x = max([elem.x for elem in chain])
    max_y = max([elem.y for elem in chain])

    return max_x - min_x, max_y - min_y


# Data structure containing the protein chain and a lattice bidirectional lookup structure
# This way super fast (neighbour) lookups can be achieved in O(1)
class ProteinLattice:

    # Initializes a new Lattice based on the given chain
    def __init__(self, chain: List[Monomer]):
        # The protein chain as a list
        self.chain: List[Monomer] = chain[:]
        self.undo_set: List[MonomerMoveRecord] = []
        # The lattice, used for fast lookups!
        self.__calculate_lattice()

    # Computes the internal lattice structure
    # We compute the grid size such that it can
    # contain the protein at it's biggest size.
    # This avoids costly re-allocations at runtime (in exchange for some memory)
    def __calculate_lattice(self):
        self.__lattice = {}

        # Set values with offset so 0,0 is in the middle of the lattice.
        for i in range(0, len(self.chain)):
            monomer = self.chain[i]
            self.__lattice[(monomer.x, monomer.y)] = MonomerRecord(MonomerRecordValue(int(monomer.kind)), i)

    # Returns an idx,value pair for a given position. idx = -1 if no monomer is present.
    def get_by_coordinate(self, x: int, y: int) -> (int, MonomerRecordValue):
        val = self.__lattice.get((x, y))
        if val is not None:
            return val.index, val.value
        return -1, MonomerRecordValue.NONE

    # Returns the monomer at position x,y
    def get_monomer(self, x: int, y: int) -> Optional[Monomer]:
        value = self.__lattice.get((x, y))
        if value is not None and value.index != -1:
            return self.chain[value.index]
        return None

    # Returns whether there is a monomer at x,y
    def has_monomer(self, x: int, y: int) -> bool:
        value = self.__lattice.get((x, y))
        if value is not None:
            return value.index != -1
        return False

    # Undoes the latest move(s) in the chain.
    def undo_last_change(self):
        # Remove all old items
        for record in self.undo_set:
            del self.__lattice[record.new]

        # Update the lattice
        for record in self.undo_set:
            monomer = self.chain[record.index]
            self.__lattice[record.old] = MonomerRecord(MonomerRecordValue(monomer.kind), record.index)
            monomer.x = record.old[0]
            monomer.y = record.old[1]

        # Clear undo set
        self.undo_set = []

    # Move monomer to different position, replacing whatever was at x,y.
    # Assumes x,y is empty!
    # Assumes monomer.x and monomer.y are set in the lattice!
    # DO NOT USE when multiple monomers need to be moved.
    # Erases and replaces the undo stack!
    def move_monomer(self, idx: int, x: int, y: int):
        monomer = self.chain[idx]

        self.undo_set = []
        self.undo_set.append(MonomerMoveRecord(idx, (monomer.x, monomer.y), (x, y)))

        del self.__lattice[(monomer.x, monomer.y)]
        self.__lattice[(x, y)] = MonomerRecord(MonomerRecordValue(monomer.kind), idx)

        self.chain[idx].x = x
        self.chain[idx].y = y

    # Debug function checking for internal inconsistencies in the lattice structure.
    def consistency_check(self):
        for idx in range(0, len(self.chain)):
            monomer = self.chain[idx]
            lat = self.__lattice[(monomer.x, monomer.y)]
            assert lat.index == idx
        assert len(self.chain) == len(self.__lattice)

    # Moves multiple monomers at once.
    # Tuple is: (index_in_chain, (x, y))
    # Erases and replaces the undo stack!
    def move_monomers(self, new_positions: List[Tuple[int, Tuple[int, int]]]):
        # Set up the undo set
        self.undo_set = []
        for idx, (x, y) in new_positions:
            monomer = self.chain[idx]
            self.undo_set.append(MonomerMoveRecord(idx, (monomer.x, monomer.y), (x, y)))

        # Remove all old items
        for idx, _ in new_positions:
            monomer = self.chain[idx]
            ox, oy = monomer.x, monomer.y
            del self.__lattice[(ox, oy)]

        # Update the lattice
        for idx, (nx, ny) in new_positions:
            monomer = self.chain[idx]
            self.__lattice[(nx, ny)] = MonomerRecord(MonomerRecordValue(monomer.kind), idx)
            monomer.x = nx
            monomer.y = ny

    # Returns the direct neighbouring Monomers around (x,y), if any
    def get_neighbours(self, x: int, y: int) -> List[Monomer]:
        neighbours = []

        for x, y in [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]:
            (idx, value) = self.get_by_coordinate(x, y)
            if idx != -1:
                neighbours.append(self.chain[idx])

        return neighbours

    # Computes the center coordinate of the protein
    def compute_center_point(self) -> Tuple[float, float]:
        min_x = min(self.chain, key=lambda e: e.x).x
        max_x = max(self.chain, key=lambda e: e.x).x
        min_y = min(self.chain, key=lambda e: e.y).y
        max_y = max(self.chain, key=lambda e: e.y).y

        dx: float = max_x - min_x
        dy: float = max_y - min_y
        return min_x + (dx / 2), min_y + (dy / 2)

    # Computes the radius of gyration
    def compute_gyration_radius(self) -> float:
        rc = self.compute_center_point()
        sum_of_squares = 0
        for monomer in self.chain:
            rk = monomer.x, monomer.y
            # Compute the rk - rc
            neg_x, neg_y = rk[0] - rc[0], rk[1] - rc[1]
            # Compute dot product with itself
            res = (neg_x ** 2) + (neg_y ** 2)
            sum_of_squares += res
        # Compute mean
        mean_value = sum_of_squares / len(self.chain)
        # Normalize
        mean_value *= 1 / len(self.chain)
        # Take the root and return
        return math.sqrt(mean_value)


# Represents collected samples from a MMC simulation
class MMCSamples:
    def __init__(self, energy: [float], gyration_radius: [float]):
        self.energy: [float] = energy
        self.gyration_radius: [float] = gyration_radius


# Represents collected samples from Annealing MMC Simulation
class AMMCSamples(MMCSamples):
    def __init__(self, energy: [float], gyration_radius: [float], heat_capacity: float, temperature: float):
        super().__init__(energy, gyration_radius)
        self.heat_capacity: float = heat_capacity
        self.temperature: [float] = temperature
