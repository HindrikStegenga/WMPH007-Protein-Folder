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

    def get_monomer(self, x: int, y: int) -> Optional[Monomer]:
        value = self.__lattice.get((x, y))
        if value is not None and value.index != -1:
            return self.chain[value.index]
        return None

    def has_monomer(self, x: int, y: int) -> bool:
        value = self.__lattice.get((x, y))
        if value is not None:
            return value.index != -1
        return False

    # Move monomer to different position, replacing whatever was at x,y.
    # Assumes x,y is empty!
    def move_monomer(self, idx: int, x: int, y: int):
        monomer = self.chain[idx]

        self.__lattice[(monomer.x, monomer.y)] = MonomerRecord(MonomerRecordValue(MonomerRecordValue.NONE), -1)
        self.__lattice[(x, y)] = MonomerRecord(MonomerRecordValue(monomer.kind), idx)

        self.chain[idx].x = x
        self.chain[idx].y = y

    # Returns the direct neighbouring Monomers around (x,y), if any
    def get_neighbours(self, x: int, y: int) -> List[Monomer]:
        neighbours = []

        for x, y in [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]:
            (idx, value) = self.get_by_coordinate(x, y)
            if idx != -1:
                neighbours.append(self.chain[idx])

        return neighbours

    # # Plots the lattice as a matrix
    # def plot(self):
    #     nparr = np.empty([len(self.__lattice), len(self.__lattice)])
    #     for x in range(0, len(self.__lattice)):
    #         for y in range(0, len(self.__lattice)):
    #             nparr[x, y] = self.__lattice[x][y].index
    #
    #     plt.matshow(nparr, cmap=plt.cm.Blues)
    #     for x in range(0, len(self.__lattice)):
    #         for y in range(0, len(self.__lattice)):
    #             c = nparr[x, y]
    #             #plt.text(x, y, str(c), va='center', ha='center')
    #     plt.show()
