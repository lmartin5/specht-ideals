from partitions import *
import random

# run file on https://trinket.io/embed/python3

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
    I = LowerOrderIdeal(9, [Partition([3, 3, 3]), Partition([4, 2, 1, 1, 1])])
    print(I)
    print(I.smaller_ideal(3))
    print(I.smaller_ideal(2))
    print(I.hilbert_series())

def get_random_partition(n):
    size = 0
    parts = []
    while size < n:
        parts.append(random.randint(1, n - size))
        size += parts[-1]
    return Partition(parts)

def add_rational_forms(pair1, pair2):
    """
    pair = (h_expr_or_poly, d) where h has integer coeffs and d is a nonnegative int.
    Returns (num_poly, D) for num_poly / (1 - t)^D in lowest possible power.
    """
    return add_many_rational_forms([pair1, pair2])

def add_many_rational_forms(terms):
    """
    terms: list of (h_expr_or_poly, d), with integer-coefficient h's.
    Returns (num_poly, D) with num_poly in ZZ[t] and minimal D after canceling (1 - t)^k.
    """
    # Normalize inputs
    polys = [Poly(h, t, domain='ZZ') for (h, _) in terms]
    ds    = [int(d) for (_, d) in terms]
    D = max(ds) if ds else 0

    # Accumulate numerator over common denominator (1 - t)^D
    num = Poly(0, t, domain='ZZ')
    for h, d in zip(polys, ds):
        k = D - d
        num += h if k == 0 else h * (ONE_MINUS_T ** k)

    # Cancel common powers of (1 - t)
    # Equivalent test: num(1) == 0  <=> divisible by (t - 1) <=> divisible by (1 - t)
    # Repeat until it no longer divides or we've exhausted the denominator.
    while D > 0 and num.eval(1) == 0:
        # exact quotient (no remainder) division by (1 - t)
        num = num.exquo(ONE_MINUS_T)
        D -= 1

    return num, D     

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

    # p = Partition([5, 5, 1, 1, 1])
    # p1 = Partition([5, 4, 2, 1, 1])
    # # P = Partitions(7)
    # # for p in P.partitions:
    # p_ser = p.hilbert_series_ii()
    # p1_ser = p1.hilbert_series_ii()

    # print(add_rational_forms((-1*p_ser[0], p_ser[1]), p1_ser))

    P = Partitions(7)
    P.show_hasse_diagram()

    for k in range(5):
        p = Partition([1 + k, 1, 1, 1, 1, 1, 1, 1])
        print(p.hilbert_series_ii())

__main__()