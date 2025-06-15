# run file on https://trinket.io/embed/python3

import networkx as nx
import matplotlib.pyplot as plt
    
class Partition:
    def __init__(self, parts):
        self.n = sum(parts)
        self.parts = parts
        self.parts.sort(reverse=True)
        self.tparts = tuple(parts)
        self.covers = [] # partitions that this one covers
        self.covered_by = [] # partitions that cover self
        self.poset_rank = 0
        
    def len(self):
        return len(self.parts)
    
    def sum(self):
        return sum(self.parts)
        
    def corner_set(self):
        cs = []
        for i in range(self.len() - 1):
            if self.parts[i + 1] < self.parts[i]:
                cs.append([i + 1, self.parts[i]])
        cs.append([self.len(), self.parts[-1]])
        return cs
        
    def dominates(self, other):
        max_len = max(len(self.parts), len(other.parts))
        s_sum = 0
        o_sum = 0
        for i in range(max_len):
            s_sum += self.parts[i] if i < len(self.parts) else 0
            o_sum += other.parts[i] if i < len(other.parts) else 0
            if s_sum < o_sum:
                return False
        return True
        
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

class Partitions:
    def __init__(self, n):
        self.n = n
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
        
    def get_partition(self, p):
        if not isinstance(p, Partition):
            p = Partition(p)
        for p1 in self.partitions:
            if p1 == p:
                return p1
        return None
        
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
