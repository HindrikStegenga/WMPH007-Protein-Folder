from enum import IntEnum
from typing import *


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


# Data structure containing the protein chain and a lattice structure
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
        grid_size = 2 * len(self.chain) + 1
        # Initialize lattice with zero values
        self.lattice = [[MonomerRecord(MonomerRecordValue.NONE, -1) for _ in range(grid_size)]
                        for _ in range(grid_size)]

        for i in range(0, len(self.chain)):
            monomer = self.chain[i]
            offset = (len(self.chain) + 1)
            self.lattice[monomer.x + offset][monomer.y + offset].index = i
            self.lattice[monomer.x + offset][monomer.y + offset].value = MonomerRecordValue(int(monomer.kind))

    def get_by_coordinate(self, x: int, y: int) -> (int, MonomerRecordValue):
        offset = int((len(self.lattice) + 1) / 2)
        val = self.lattice[x + offset][y + offset]

        return val.index, val.value
