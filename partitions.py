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
        self.parts = list(parts)
        self.parts.sort(reverse=True)
        self.parts = [part for part in self.parts if part > 0]
        self.tparts = tuple(self.parts)

    def __len__(self):
        return len(self.parts)
    
    def __iter__(self):
        return iter(self.parts)

    # allows for accessing parts in an easier manner, i.e. p[1] will give first part of p
    # WARNING: p[len(p) + 1] will return 0, not an error
    def __getitem__(self, k):
        if isinstance(k, int):
            if k < 1:
                raise IndexError("Partition index out of range")
            if k > len(self):
                return 0
            else:
                return self.parts[k - 1]
        else:
            raise TypeError("Partition indices must be integers")
    
    # this creates a copy of the partition with part updated (done to prevent breaking hashing behavior)
    def update_part(self, k, value):
        if k < 1:
            raise IndexError("Partition index out of range")
        if value < 0:
            raise ValueError("Partition parts must be non-negative")

        new_partition = self.copy()

        # deals with case when part is 0
        if value == 0:
            if self[k] == 0:
                return new_partition
            elif k == 1 and self == Partition([1]):
                raise ValueError("Empty partition is not allowed")
            else:
                del new_partition.parts[k - 1]
                new_partition.parts.sort(reverse=True)
                new_partition.tparts = tuple(new_partition.parts)
                return new_partition

        if k > len(self):
            new_partition.parts.append(value)
        else:
            new_partition.parts[k - 1] = value
        new_partition.parts.sort(reverse=True)
        new_partition.tparts = tuple(new_partition.parts)
        return new_partition

    def copy(self):
        return Partition(self.parts.copy())

    # creates a new partition where k-th part is increased by 1
    def add_to_part(self, k):
        return self.update_part(k, self[k] + 1)
    
    # creates a new partition where k-th part is decreased by 1
    def remove_from_part(self, k):
        if k > len(self):
            raise ValueError("Cannot remove from a zero part")
        return self.update_part(k, self[k] - 1)
    
    # returns the list of rows that have a corner
    # if self = [4, 2, 2, 2], corner_set returns [1, 4]
    def corner_set(self):
        cs = []
        for i in range(1, len(self) + 1):
            if self[i] > self[i + 1]:
                cs.append(i)
        return cs

    def compress(self, k):
        if k > len(self) or k < 1:
            raise ValueError("Can only compress on non-zero row")
        if self[k] < 2:
            raise ValueError("Cannot compress on part of size 1")
            
        compressed_partition = self.copy()
        m = max(i for i in range(1, len(self) + 1) if self[i] == self[k])
        blocks_to_move = m - k

        # we have to update the k + 1 part every time because this is where the part with size self[k] will be after reordering
        for i in range(k + 1, m + 1):
            compressed_partition = compressed_partition.update_part(k + 1, self[k] - 1)

        i = m + 1
        while blocks_to_move > 0:
            compressed_partition = compressed_partition.update_part(i, min(self[k] - 1, self[i] + blocks_to_move))
            blocks_to_move -= self[k] - 1 - self[i]
            i += 1

        return compressed_partition
        
    def dominates(self, other):
        s_sum = 0
        o_sum = 0
        index = 1
        while self[index] > 0 or other[index] > 0:
            s_sum += self[index]
            o_sum += other[index]
            index += 1
            if o_sum < s_sum:
                return False
        return True
    
    def meet(self, other):
        max_len = max(len(self), len(other))
        s_sum = 0
        o_sum = 0
        meet_parts = []
        m_sum = 0
        for i in range(1, max_len + 1):
            s_sum += self[i]
            o_sum += other[i]
            m_part = min(s_sum, o_sum) - m_sum
            m_sum += m_part

            if m_part == 0:
                return Partition(meet_parts)
            else:
                meet_parts.append(m_part)

    def less_than_ideal(self):
        return LowerOrderIdeal(sum(self), [self])
    
    def strictly_less_than_ideal(self):
        return LowerOrderIdeal(sum(self), self.get_covered_partitions())
        
    def get_covered_partitions(self):
        covered_partitions = []
        for i in self.corner_set():
            j = self[i]
            if i == len(self):
                if j > 1:
                    new_partition = self.update_part(i, j - 1)
                    new_partition = new_partition.update_part(i + 1, 1)
                    covered_partitions.append(new_partition)
            elif self[i + 1] < j - 1:
                new_partition = self.update_part(i, j - 1)
                new_partition = new_partition.update_part(i + 1, self[i + 1] + 1)
                covered_partitions.append(new_partition)
            elif j == 2:
                new_partition = self.update_part(i, j - 1)
                new_partition = new_partition.update_part(len(new_partition) + 1, 1)
                covered_partitions.append(new_partition)
            else:
                for k in range(i, len(self) + 1):
                    if self[k] == j - 2:
                        new_partition = self.update_part(i, j - 1)
                        new_partition = new_partition.update_part(k, self[k] + 1)
                        covered_partitions.append(new_partition)
                        break
        return covered_partitions
    
    def get_chipping_sequences(self):
        if self[1] == 1:
            return [[self]]
        
        sequences = []
        for i in self.corner_set():
            new_sequences = self.remove_from_part(i).get_chipping_sequences()
            for seq in new_sequences:
                sequences.append([self] + seq)
        
        first_corner = min(self.corner_set())
        for i in range(1, len(self) + 1):
            if self[i] == 1:
                break
            if i in self.corner_set():
                continue

            ith_compression = self.compress(i)
            new_sequences = ith_compression.remove_from_part(i).get_chipping_sequences()
            for seq in new_sequences:
                sequences.append([self] + seq)

            if i >= first_corner:
                ci = max(j for j in self.corner_set() if j <= i)
                new_sequences = ith_compression.remove_from_part(ci).get_chipping_sequences()
                for seq in new_sequences:
                    sequences.append([self] + seq)

        return sequences

    # allows for use of Partitions in sets and dicts
    def __hash__(self):
        return hash(self.tparts)
        
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
        return f"{self.tparts}"
    
### ------------------------------------------------------------------------ ###

class PartitionSet:
    def __init__(self, partitions=None):
        if partitions is None:
            partitions = []
        self.partitions = partitions

    def __len__(self):
        return len(self.partitions)
    
    def __iter__(self):
        return iter(self.partitions)
        
    def get_partition(self, p):
        if not isinstance(p, Partition):
            p = Partition(p)
        for p1 in self:
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
        return f"Set of {len(self)} partitions: {self.partitions}"
     
### ------------------------------------------------------------------------ ###

class Partitions(PartitionSet):
    def __init__(self, n):
        super().__init__()
        if n < 1:
            raise ValueError("In PartitionSet, n must be a positive integer")
        
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
                covered_partitions = p.get_covered_partitions()
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
    # with current setup, the ideals are allowed to be empty, when no generators are specified
    def __init__(self, n, generators=None):
        super().__init__()
        if n < 1:
            raise ValueError("In PartitionSet, n must be a positive integer")
        
        self.n = n
        if generators is None:
            generators = []

        self.generators = []
        self.partitions = []
        for g in generators:
            self.add_partition_to_ideal(g)
    
    # add partition and all lower partitions
    def add_partition_to_ideal(self, p):
        if sum(p) != self.n:
                raise ValueError("All partitions of the ideal must be partitions of n = " + str(self.n))
        
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

            partitions_to_add = p.get_covered_partitions()
        else:
            return False
        
        while len(partitions_to_add) > 0:
            p_to_add = partitions_to_add[0]
            success = self.add_partition(p_to_add)
            if success:
                covered_partitions = p_to_add.get_covered_partitions()
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

p = Partition([8, 5, 1, 1])
cs = p.get_chipping_sequences()
for s in cs:
    print(s)
