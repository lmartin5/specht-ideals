from partitions import *
import random

def partition_example():
    p1 = Partition([4, 2, 2, 2])
    p2 = Partition([5, 2, 1, 1, 1])
    print(p1, p2)
    print(p1 >= p2, p1 <= p2, p1 >= Partition([1] * 10))
    print(p1.corner_set(), p2.corner_set())
    print(p1.add_to_part(3), p2.add_to_part(2), p2.add_to_part(7))

# TODO: allow for empty lower order ideals, should fix my issue?
def partitions_and_ideals_example():
    P = Partitions(8)
    print(P)
    # P.show_hasse_diagram()
    print(len(P.partitions))
    I = LowerOrderIdeal(9, [Partition([3, 3, 3])])
    print(I)
    print(I.smaller_ideal(3))
    print(I.smaller_ideal(2))

def get_random_partition(n):
    size = 0
    parts = []
    while size < n:
        parts.append(random.randint(1, n - size))
        size += parts[-1]
    return Partition(parts)        

def __main__():
    # example of using partitions
    # partition_example()
    # print("**********************")
    # partitions_and_ideals_example()
    
    # while True:
    #     p = get_random_partition(15)
    #     r_row = random.randrange(0, p.len())
    #     r_part = p.parts[r_row]

    #     corner_rows = [x[0] for x in p.corner_set()]
    #     m = max((x for x in corner_rows if x <= r_row + 1), default=0)

    #     if r_part == 1 or m == 0 or m == r_row + 1:
    #         continue
    #     else:
    #         break

    # print("p:", p)
    # print("Row:", r_row + 1)
    # print("C_t:", m)
    # p_1 = Partition([x - 1 if idx == m - 1 else x for idx, x in enumerate(p.parts)])
    # gamma_parts = p.parts[0: r_row + 1]
    # gamma_parts += ([r_part - 1] * ((p.sum() - sum(gamma_parts)) // (r_part - 1)))
    # if (p.sum() - sum(gamma_parts)) % (r_part - 1) != 0:
    #     gamma_parts.append((p.sum() - sum(gamma_parts)) % (r_part - 1))

    # mt = Partition(gamma_parts).meet(p)
    # p_2 = Partition([x - 1 if idx == r_row else x for idx, x in enumerate(mt.parts)])

    # print("p1:", p_1)
    # print("p2:", p_2)
    # print("meet:", p_1.meet(p_2))

    p = Partition([9, 5, 5, 2])
    #P = Partitions(4)
    #for p in P.partitions:
    #    print(p.tparts, ":", sym.simplify(p.hilbert_series_ii()))
    # print("First:", p.hilbert_series())
    print(sym.simplify(p.hilbert_series_ii()))

__main__()