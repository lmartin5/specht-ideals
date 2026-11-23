from math import comb
from partitions import *

def f(a, b, c, k):
    return (
        k*(a+b+c-2)
        + comb(a+b+c,2)
        + (a+b)*(2*c+b)
        + (c-1)*b + 4*comb(c,2) - 3
    )

def g(a, b, c, k):
    return (
        (k+1)*(b+c-2)
        + 3*comb(b,2)
        + 5*comb(c,2)
        + 4*b*c + b + c - 4
    )

def h(a, b, c, k):
    return (
        k*(b-2)
        + 3*comb(b,2)
        + b - 3
    )

# example: build matrix for fixed b, fixed c
def matrix(b, c, a_max=9, k_max=7):
    M = []
    for a in range(1, a_max+1):
        row = []
        for k in range(0, k_max+1):
            row.append(f(a,b,c,k))
        M.append(row)
    return M

def formula(p):
    w = len(p)
    a = p.parts.count(1)

    deg = 0
    deg += comb(w,2)
    deg += (w - 2) * p[1]
    deg += (w - 1) * (sum(p[i] for i in range(2, w - a + 1)) - (w - a))
    
    return deg

# a, b, c = 0, 0, 0
# for a in range(0,1):
#     for b in range(2,7):
#         info = []
#         for k in range(0,8):
#             info.append([fun1(a,b,k), f(a,b,0,k)])
#         print("a, b, c:", a, b, c, "degrees:", info)
#     print("**********************")

# c = 5
# for b in range(4):
#     print("b:", b)
#     for row in matrix(b, c):
#         print(row)

for p in Partitions(15):
    print(p, formula(p))
    print("****************")