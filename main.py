from generation import *
from drawing import *
import math
import numpy


def boltzmann(x: float, temperature: float) -> float:
    return math.exp(- x / float(1.0) * temperature)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    samples = mmc(5.0, 25, 5000, 100, 0.5)

    plt.plot(samples)
    plt.show()
