% https://python-fiddle.com/examples/matplotlib?checkpoint=1749747241

import networkx as nx
import matplotlib.pyplot as plt

def draw_hasse(poset, cover_relations, label_func=str):
    """
    Draws the Hasse diagram of a poset.

    Parameters:
        poset: A list of elements.
        cover_relations: A list of (a, b) tuples where `a < b` and `b` covers `a`.
        label_func: A function to convert elements to displayable labels.
    """
    G = nx.DiGraph()
    G.add_nodes_from(poset)
    G.add_edges_from(cover_relations)

    # Use graphviz layout if available for better Hasse look

    pos = nx.spring_layout(G)

    labels = {node: label_func(node) for node in G.nodes()}

    plt.figure(figsize=(10, 20))
    nx.draw(G, pos, labels=labels, with_labels=True, node_size=1500,
            node_color='lightblue', font_size=10, arrows=False)

    plt.title("Hasse Diagram")
    plt.axis('off')
    plt.show()
    
class Partition:
    def __init__(self, parts):
        self.parts = parts
        self.parts.sort(reverse=True)
        self.tparts = tuple(parts)
        self.covered_by = []   # Partitions that cover this one
        self.covered = []  # Partitions that this one covers
        self.rank = 0
        
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
        
    def __hash__(self):
        return hash(self.tparts)  # allows use in sets and dicts
        
    def __eq__(self, other):
        if not isinstance(other, Partition):
            return False
        return self.tparts == other.tparts

    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return self._dominates(other)

    def __ge__(self, other):
        return other._dominates(self)

    def __lt__(self, other):
        return self <= other and self != other

    def __gt__(self, other):
        return self >= other and self != other

    def _dominates(self, other):
        max_len = max(len(self.parts), len(other.parts))
        s_cumsum = 0
        o_cumsum = 0
        for i in range(max_len):
            s_cumsum += self.parts[i] if i < len(self.parts) else 0
            o_cumsum += other.parts[i] if i < len(other.parts) else 0
            if s_cumsum < o_cumsum:
                return False
        return True

    def __repr__(self):
        return f"Partition({self.parts})"

class Partitions:
    def __init__(self, n):
        self.n = n
        self.partitions = self.generate_partitions()

    def generate_partitions(self):
        def find_covered_partitions(p):
            covered = []
            for i, j in p.corner_set():
                if i == p.len():
                    if j > 1:
                        new_parts = p.parts.copy()
                        new_parts[-1] += -1
                        new_parts.append(1)
                        covered.append(Partition(new_parts))
                elif p.parts[i] < j - 1:
                    new_parts = p.parts.copy()
                    new_parts[i - 1] += -1
                    new_parts[i] += 1
                    covered.append(Partition(new_parts))
                elif j == 2:
                    new_parts = p.parts.copy()
                    new_parts[i - 1] += -1
                    new_parts.append(1)
                    covered.append(Partition(new_parts))
                else:
                    for k in range(i, p.len()):
                        if p.parts[k] == j - 2:
                            new_parts = p.parts.copy()
                            new_parts[i - 1] += -1
                            new_parts[k] += 1
                            covered.append(Partition(new_parts))
                            break
            return covered
            
        max_partition = Partition([self.n])
        max_partition.rank = 0
        partitions = [max_partition]
        new_partitions = [max_partition]
        new_partitions_extra = []
        
        while len(new_partitions) != 0: 
            for p in new_partitions:
                covered_partitions = find_covered_partitions(p)
                p.covered = covered_partitions
                for p1 in covered_partitions:
                    add_to_list = True
                    
                    for p2 in partitions:
                        if p1.parts == p2.parts:
                            p1 = p2
                            add_to_list = False
                            break
                    
                    p1.covered_by.append(p)
                    p1.rank = max(p1.rank, p.rank + 1)
                    
                    if add_to_list: 
                        partitions.append(p1)
                        new_partitions_extra.append(p1)
                        
            new_partitions = new_partitions_extra
            new_partitions_extra = []
            
        # for p in partitions:
        #    p.rank = max([p1.rank if p1.parts != [p.sum()] else 0 for p1 in p.covered_by]) + 1
            
        return partitions
        
    def get_cover_relations(self):
        relations = []
        for p in self.partitions:
            for p1 in p.covered:
                relations.append((p1, p))
        return relations

    def __repr__(self):
        return f"Partitions({self.n}) with {len(self.partitions)} partitions"

def draw_ranked_hasse(poset, cover_relations, label_func=str, rank_func=None):
    import networkx as nx
    import matplotlib.pyplot as plt

    G = nx.DiGraph()
    G.add_nodes_from(poset)
    G.add_edges_from(cover_relations)

    # Use a manual layout to control vertical levels
    if rank_func is None:
        # Default rank: negate length of partition (more parts = lower rank)
        rank_func = lambda p: -len(p.parts)

    levels = {}
    for node in G.nodes():
        levels.setdefault(rank_func(node), []).append(node)

    pos = {}
    for y, level in sorted(levels.items(), reverse=True):
        n = len(level)
        if n == 1:
            # Single node: put exactly in the middle (x=0)
            pos[level[0]] = (0, y)
        else:
            # Multiple nodes: space evenly centered at x=0
            # For example: for n=3, x = [-1, 0, 1]
            spacing = 2  # distance between nodes
            xs = [spacing * (i - (n - 1) / 2) for i in range(n)]
            for x, node in zip(xs, level):
                pos[node] = (x, y)

    labels = {node: label_func(node) for node in G.nodes()}

    plt.figure(figsize=(10, 20))
    nx.draw(G, pos, labels=labels, with_labels=True, node_size=1500,
            node_color='lightblue', font_size=10, arrows=False)

    plt.title("Hasse Diagram (Custom Layout)")
    plt.axis('off')
    plt.show()

P = Partitions(10)
poset = P.partitions
cover_relations = P.get_cover_relations()

draw_ranked_hasse(poset, cover_relations, lambda p: str(p.parts), lambda p: -p.rank)
