'''
File contains:
    class Partition
    class PartitionSet
    class Partitions
    class LowerOrderIdeal
'''
import networkx as nx
import matplotlib.pyplot as plt

class Partition:
    # parts should be a list of positive numbers, i.e. [3, 2, 2, 1]
    def __init__(self, parts):
        self.n = sum(parts)
        self.parts = parts
        self.parts.sort(reverse=True)
        self.tparts = tuple(parts)
        
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
    def __init__(self, partitions=None):
        if partitions is None:
            partitions = []
        self.partitions = partitions
        
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
        return f"Set of {len(self.partitions)} partitions: {self.partitions}"
     
### ------------------------------------------------------------------------ ###

class Partitions(PartitionSet):
    def __init__(self, n):
        super().__init__()
        self.n = n
        self.covers = {} # dict holding covering relations
        self.covered_by = {} # dict holding covering relations
        self.ranks = {} # dict holding poset ranks to help generate the Hasse diagram
        self.generate_partitions()

    def generate_partitions(self):
        max_partition = Partition([self.n])
        self.add_partition(max_partition)
        self.covered_by[max_partition] = []
        self.ranks[max_partition] = 0

        new_partitions = [max_partition]
        new_partitions_extra = []
        
        while len(new_partitions) != 0:
            for p in new_partitions:
                covered_partitions = p.find_covered_partitions()
                self.covers[p] = covered_partitions
                for p1 in covered_partitions:
                    # add partition if not already added
                    added_to_list = self.add_partition(p1)
                    if added_to_list:
                        self.ranks[p1] = 0
                        self.covered_by[p1] = []
                        new_partitions_extra.append(p1)
                        
                    self.covered_by[p1].append(p)
                    self.ranks[p1] = max(self.ranks[p1], self.ranks[p] + 1)
                        
            new_partitions = new_partitions_extra
            new_partitions_extra = []
        
    def get_cover_relations(self):
        relations = []
        for p in self.partitions:
            for p1 in self.covers[p]:
                relations.append((p1, p))
        return relations

    def show_hasse_diagram(self, width=10, height=20):
        G = nx.DiGraph()
        G.add_nodes_from(self.partitions)
        G.add_edges_from(self.get_cover_relations())
    
        # using a manual layout to control vertical levels
        levels = {}
        for node in G.nodes():
            levels.setdefault(-1*self.ranks[node], []).append(node)
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
    # generators should be a list of partitions of n
    # with current setup, the ideals are allowed to be empty
    def __init__(self, n, generators=None):
        super().__init__()
        self.n = n
        if generators is None:
            generators = []
        
        for g in generators:
            if g.sum() != n:
                raise ValueError("All generators of the ideal must be partitions of n = " + str(n))

        self.generators = []
        self.partitions = []
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

            partitions_to_add = p.find_covered_partitions()
        else:
            return False
        
        while len(partitions_to_add) > 0:
            p_to_add = partitions_to_add[0]
            success = self.add_partition(p_to_add)
            if success:
                covered_partitions = p_to_add.find_covered_partitions()
                partitions_to_add.extend(covered_partitions)
            partitions_to_add.pop(0)
        return True

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
    
### ------------------------------------------------------------------------ ###