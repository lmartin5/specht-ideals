# run file on https://trinket.io/embed/python3
import networkx as nx
import matplotlib.pyplot as plt
import sympy as sym
from sympy import symbols, Poly
import math as math

t = symbols('t')
ONE_MINUS_T = Poly(1 - t, t, domain='ZZ')

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

class Partition:
    # parts should be a list of positive numbers, i.e. [3, 2, 2, 1]
    def __init__(self, parts):
        self.n = sum(parts)
        self.parts = parts
        self.parts.sort(reverse=True)
        self.tparts = tuple(parts)
        self.covers = [] # partitions that self covers
        self.covered_by = [] # partitions that cover self
        self.poset_rank = 0 # keeps track of where partition should be placed in Hasse diagram
        
    def len(self):
        return len(self.parts)
    
    def sum(self):
        return sum(self.parts)
    
    # if self = [4, 2, 2, 2], corner_set returns [[1, 4], [4, 2]]
    def corner_set(self):
        cs = []
        for i in range(self.len() - 1):
            if self.parts[i + 1] < self.parts[i]:
                cs.append([i + 1, self.parts[i]])
        cs.append([self.len(), self.parts[-1]])
        return cs
    
    # creates a new partition where k-th part is increased by 1
    def add_to_part(self, k):
        new_parts = self.parts.copy()
        if k >= self.len() + 1:
            new_parts.append(1)
        else:
            new_parts[k - 1] += 1
        return Partition(new_parts)
    
    def remove_from_part(self, k):
        new_parts = self.parts.copy()
        if k > self.len():
            k = self.len()

        if new_parts[k - 1] == 1:
            new_parts.remove(1)
        else:
            new_parts[k - 1] += -1
        return Partition(new_parts)      
        
    def dominates(self, other):
        min_len = min(len(self.parts), len(other.parts))
        s_sum = 0
        o_sum = 0
        for i in range(min_len):
            s_sum += self.parts[i]
            o_sum += other.parts[i]
            if o_sum < s_sum:
                return False
        return True
    
    def meet(self, other):
        max_len = max(len(self.parts), len(other.parts))
        s_sum = 0
        o_sum = 0
        meet_parts = []
        m_sum = 0
        for i in range(max_len + 1):
            s_sum += 0 if i >= len(self.parts) else self.parts[i]
            o_sum += 0 if i >= len(other.parts) else other.parts[i]
            m_part = min(s_sum, o_sum) - m_sum
            m_sum += m_part

            if m_part == 0:
                return Partition(meet_parts)
            else:
                meet_parts.append(m_part)

        
    def find_covered_partitions(self):
        p = self
        covered_partitions = []
        for i, j in p.corner_set():
            if i == p.len():
                if j > 1:
                    new_parts = p.parts.copy()
                    new_parts[-1] += -1
                    new_parts.append(1)
                    covered_partitions.append(Partition(new_parts))
            elif p.parts[i] < j - 1:
                new_parts = p.parts.copy()
                new_parts[i - 1] += -1
                new_parts[i] += 1
                covered_partitions.append(Partition(new_parts))
            elif j == 2:
                new_parts = p.parts.copy()
                new_parts[i - 1] += -1
                new_parts.append(1)
                covered_partitions.append(Partition(new_parts))
            else:
                for k in range(i, p.len()):
                    if p.parts[k] == j - 2:
                        new_parts = p.parts.copy()
                        new_parts[i - 1] += -1
                        new_parts[k] += 1
                        covered_partitions.append(Partition(new_parts))
                        break
        return covered_partitions
    
    def hilbert_series_ii(self):
        def hilbert_series_finder(p, series):
            if p.parts[0] == 1:
                return (1 - t**(math.comb(p.n, 2))), p.n
            if p.len() == 1:
                return 0, 0
            
            # recursion
            new_terms = []
            for i in range(p.len() - 1):
                Li_HS = Li_hilbert_series_finder(p, i + 1, series)
                new_terms.append((t**(i) * Li_HS[0], Li_HS[1]))
            Li_HS = Li_hilbert_series_finder(p, p.len(), series)
            new_terms.append((t**(p.len() - 1) * Li_HS[0], Li_HS[1] + 1))

            return add_many_rational_forms(new_terms)

        def Li_hilbert_series_finder(p, j, series):
            part = p.parts[j - 1]
            corner_rows = [x[0] for x in p.corner_set()]
            m = max((x for x in corner_rows if x <= j), default=0)

            if part == 1:
                p1 = p.remove_from_part(m)
                return series[p1.tparts]

            gamma_parts = p.parts[0: j]
            gamma_parts += ([part - 1] * ((p.sum() - sum(gamma_parts)) // (part - 1)))
            if (p.sum() - sum(gamma_parts)) % (part - 1) != 0:
                gamma_parts.append((p.sum() - sum(gamma_parts)) % (part - 1))

            mt = Partition(gamma_parts).meet(p)
            p2 = mt.remove_from_part(j)

            if m == 0:
                return series[p2.tparts]
            
            p1 = p.remove_from_part(m)
            return add_many_rational_forms([series[p1.tparts], series[p2.tparts], (-1*series[p1.meet(p2).tparts][0], series[p1.meet(p2).tparts][1])])

        series = {}

        for i in range(2, self.n):
            print(i)
            P = Partitions(i)
            for p in P.partitions:
                series[p.tparts] = hilbert_series_finder(p, series)
        
        return hilbert_series_finder(self, series)
    
    def hilbert_series(self):
        return sym.simplify(self.hilbert_recursion())
    
    def hilbert_recursion(self):
        # base cases
        if self.parts[0] == 1:
            return (1 - t**(math.comb(self.n, 2))) / (1 - t)**self.n
        if self.len() == 1:
            return 0

        # recursion
        hs = 0
        for i in range(self.len() - 1):
            hs += t**(i) * self.get_Li_hilbert_series(i + 1)
        hs += (t**(self.len() - 1) / (1 - t)) * self.get_Li_hilbert_series(self.len())
        return hs
    
    def get_Li_hilbert_series(self, i):
        part = self.parts[i - 1]
        corner_rows = [x[0] for x in self.corner_set()]
        m = max((x for x in corner_rows if x <= i), default=0)

        if part == 1:
            p1 = self.remove_from_part(m)
            return p1.hilbert_recursion()

        gamma_parts = self.parts[0: i]
        gamma_parts += ([part - 1] * ((self.sum() - sum(gamma_parts)) // (part - 1)))
        if (self.sum() - sum(gamma_parts)) % (part - 1) != 0:
            gamma_parts.append((self.sum() - sum(gamma_parts)) % (part - 1))

        mt = Partition(gamma_parts).meet(self)
        p2 = mt.remove_from_part(i)

        if m == 0:
            return p2.hilbert_recursion()
        
        p1 = self.remove_from_part(m)
        return (p1.hilbert_recursion() + p2.hilbert_recursion() - p1.meet(p2).hilbert_recursion())

        
    def __hash__(self):
        return hash(self.tparts) # allows for use in sets and dicts
        
    def __eq__(self, other):
        if not isinstance(other, Partition):
            return False
        return self.tparts == other.tparts

    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return self.dominates(other)

    def __ge__(self, other):
        return other.dominates(self)

    def __lt__(self, other):
        return self <= other and self != other

    def __gt__(self, other):
        return self >= other and self != other

    def __repr__(self):
        return f"Partition({self.parts})"

### ------------------------------------------------------------------------ ###

class PartitionSet:
    def __init__(self, n):
        self.n = n
        self.partitions = []
        
    def get_partition(self, p):
        if not isinstance(p, Partition):
            p = Partition(p)
        for p1 in self.partitions:
            if p1 == p:
                return p1
        return False
        
    def add_partition(self, p):
        if not self.get_partition(p):
            self.partitions.append(p)
            return True
        else:
            return False
            
    def remove_partition(self, p):
        p1 = self.get_partition(p)
        if not p1:
            return False
        else:
            self.partitions.remove(p1)
            return True

    def __repr__(self):
        return f"Set of {len(self.partitions)} partitions of {self.n}"
     
### ------------------------------------------------------------------------ ###

class Partitions(PartitionSet):
    def __init__(self, n):
        super().__init__(n)
        self.partitions = self.generate_partitions()

    def generate_partitions(self):
        max_partition = Partition([self.n])
        partitions = [max_partition]
        new_partitions = [max_partition]
        new_partitions_extra = []
        
        while len(new_partitions) != 0: 
            for p in new_partitions:
                covered_partitions = p.find_covered_partitions()
                p.covers = covered_partitions
                for p1 in covered_partitions:
                    # check if partition is already added to list of partitions
                    add_to_list = True
                    for p2 in partitions:
                        if p1.parts == p2.parts:
                            p1 = p2
                            add_to_list = False
                            break
                    if add_to_list: 
                        partitions.append(p1)
                        new_partitions_extra.append(p1)
                        
                    p1.covered_by.append(p)
                    p1.poset_rank = max(p1.poset_rank, p.poset_rank + 1)
                        
            new_partitions = new_partitions_extra
            new_partitions_extra = []
            
        return partitions
        
    def get_cover_relations(self):
        relations = []
        for p in self.partitions:
            for p1 in p.covers:
                relations.append((p1, p))
        return relations

    def show_hasse_diagram(self, width=10, height=20):
        G = nx.DiGraph()
        G.add_nodes_from(self.partitions)
        G.add_edges_from(self.get_cover_relations())
    
        # using a manual layout to control vertical levels
        levels = {}
        for node in G.nodes():
            levels.setdefault(-1*node.poset_rank, []).append(node)
        pos = {}
        for y, level in sorted(levels.items(), reverse=True):
            l = len(level)
            if l == 1:
                # Single node: put exactly in the middle (x=0)
                pos[level[0]] = (0, y)
            else:
                # Multiple nodes: space evenly centered at x=0
                # For example: for n=3, x = [-1, 0, 1]
                spacing = 2  # distance between nodes
                xs = [spacing * (i - (l - 1) / 2) for i in range(l)]
                for x, node in zip(xs, level):
                    pos[node] = (x, y)

        labels = {node: str(node.parts) for node in G.nodes()}
        plt.figure(figsize=(width, height))
        nx.draw(G, pos, labels=labels, with_labels=True, node_size=1500,
                node_color='lightblue', font_size=10, arrows=False)
        plt.title(f'Hasse Diagram for $P_{{{self.n}}}$')
        plt.axis('off')
        plt.show()

    def __repr__(self):
        return f"Partitions({self.n}) with {len(self.partitions)} partitions"
        
### ------------------------------------------------------------------------ ###

class LowerOrderIdeal(PartitionSet):
    # g should be a list of partitions of n
    def __init__(self, n, generators=[]):
          super().__init__(n)
          self.generators = []
          self.partitions = [Partition([1] * (self.n))]
          for g in generators:
            self.add_partition_to_ideal(g)
    
    # add partition and all lower partitions
    def add_partition_to_ideal(self, p):
        success = self.add_partition(p)
        if success:
            self.generators.append(p)
            # remove unnecessary generators
            unnec_gens = []
            for g in self.generators:
                if g < p:
                    unnec_gens.append(g)
            for g in unnec_gens:
                self.generators.remove(g)
            
            covered_partitions = p.covers.copy()
            if len(covered_partitions) == 0:
                covered_partitions = p.find_covered_partitions()
            partitions_to_add = covered_partitions
        else:
            return False
        
        while len(partitions_to_add) > 0:
            p_to_add = partitions_to_add[0]
            success = self.add_partition(p_to_add)
            if success:
                covered_partitions = p_to_add.covers.copy()
                if len(covered_partitions) == 0:
                    covered_partitions = p_to_add.find_covered_partitions()
                partitions_to_add.extend(covered_partitions)
            partitions_to_add.pop(0)
        return True
        
    def hilbert_series(self):
        return sym.simplify(self.hilbert_recursion())
        
    def hilbert_recursion(self):
        # base cases
        if self.n == 1 or (len(self.generators) == 1 and self.generators[0] == Partition([self.n])):
            return 0
            
        if len(self.generators) == 0 or (len(self.generators) == 1 and self.generators[0] == Partition([1] * self.n)):
            return (1 - t**(math.comb(self.n, 2))) / (1 - t)**self.n
            
        if len(self.generators) == 1 and self.generators[0] == Partition([self.n - 1] + [1]):
            return 1 / (1 - t)
            
        if len(self.generators) == 1 and self.generators[0].len() == 2:
            n = self.n
            m = self.generators[0].parts[1]
            numer = 1
            for i in range(1, m):
                numer += math.comb(n - m + i - 1, i) * t**i
            numer += math.comb(n - 1, m - 2) * t**m
            return numer / (1 - t)**m
            
        if len(self.generators) == 1 and self.generators[0].len() == 3 and self.generators[0].parts[0] == self.generators[0].parts[1] and self.generators[0].parts[2] == 1:
            n = self.n
            m = self.generators[0].parts[0]
            numer = 1
            for i in range(1, m + 2):
                numer += math.comb(m + i - 1, i) * t**i
            return numer / (1 - t)**m
        
        # recursion
        hs = 0
        max_len = max([g.len() for g in self.generators])
        for i in range(max_len - 1):
            hs += t**(i) * self.create_ideal_for_smaller_n(i + 1).hilbert_recursion()
        hs += (t**(max_len - 1) / (1 - t)) * self.create_ideal_for_smaller_n(max_len).hilbert_recursion()
        return hs
    
    # returns the lower order ideal of P_{n - 1} of all partitions that are in I after adding 1 to the k-th part
    def smaller_ideal(self, k):
        P = Partitions(self.n - 1)
        I = LowerOrderIdeal(self.n - 1)
        
        for p in P.partitions:
            for g in self.generators:
                if p.add_to_part(k) <= g:
                    I.add_partition_to_ideal(p)
                    break
        
        return I
        
    def __repr__(self):
        return f"Lower order ideal of {len(self.partitions)} partitions of {self.n}\nGenerated by {self.generators}"
