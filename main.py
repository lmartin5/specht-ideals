# run file on https://trinket.io/embed/python3
from partitions import *
from specht_ideals import specht_hilbert_series, ideal_specht_hilbert_series

def __main__():
    # partition examples
    p = Partition([4, 2, 1])
    q = Partition([3, 3, 1])
    print(p)
    print(p.corner_set())
    print(p >= q) # checked in the dominance order
    print(p.meet(q)) # meet in the dominance order

    # full poset examples
    P = Partitions(5)
    print(P)
    P.show_hasse_diagram()

    # lower order ideal examples
    g = Partition([3, 2])
    I = LowerOrderIdeal(5, generators=[g])
    print(I)
    print(I.partitions)
    J = I.smaller_ideal(1)
    print(J)
    print(J.partitions)

    # hilbert series examples
    print(specht_hilbert_series(p))
    print(ideal_specht_hilbert_series(I))

__main__()