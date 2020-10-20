from generation import *
from drawing import *


# Checks whether the random walk creates proteines correctly.
# Only checks the hydrophobicity fraction!
def check_random_walk():
    ratios = []
    for i in range(0, 25):
        protein = generate_protein(25, 0.5)
        hcount = get_kind_count(protein, MonomerKind.H)
        pcount = get_kind_count(protein, MonomerKind.P)
        ratio = float(hcount) / 25.0
        ratios.append(ratio)
        print('{} {}/{} - {}'.format(ratio,
                                     hcount,
                                     pcount,
                                     get_chain_composition_string(protein)))

    print('Average hydrophobicity: {:.2f}'.format(sum(ratios) / len(ratios)))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    check_random_walk()
    samples = mmc(5.0, 25, 5000, 100, 0.5)
    plt.plot(samples)
    plt.show()
