from matplotlib import pyplot as plt
import numpy as np
from enum import IntEnum
from typing import *


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

    # Move monomer to different position, replacing whatever was at x,y.
    # Assumes x,y is empty!
    # Assumes monomer.x and monomer.y are set in the lattice!
    # DO NOT USE when multiple monomers need to be moved.
    def move_monomer(self, idx: int, x: int, y: int):
        monomer = self.chain[idx]

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
    def move_monomers(self, new_positions: List[Tuple[int, Tuple[int, int]] ]):
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
