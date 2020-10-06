from enum import IntEnum


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