import os.path
import os
from pyvis.network import Network
import networkx as nx


class Graph:
    def __init__(self):
        # Initialize an empty dictionary to store the graph
        self.graph = nx.DiGraph()
        self.label2idx = {}

    def add_vertex(self, vertex):
        # Add a vertex to the graph
        if vertex not in self.label2idx.keys():
            current_idx = len(self.label2idx.keys())
            self.graph.add_node(current_idx, label=vertex)
            self.label2idx[vertex] = current_idx

    def add_edge(self, vertex1, vertex2):
        # Add an edge between two vertices
        v1_idx = self.label2idx[vertex1]
        v2_idx = self.label2idx[vertex2]
        self.graph.add_edge(v1_idx, v2_idx)

    def calculate_degree(self, vertex):
        vertex_idx = self.label2idx[vertex]
        return self.graph.degree[vertex_idx]

    def get_vertices(self):
        # Get all vertices in the graph
        return self.graph.nodes

    def get_edges(self):
        # Get all edges in the graph
        return self.graph.edges

    def plot(self, width=500, height=500):
        net = Network('1000px', '1000px', directed=True)
        net.from_nx(self.graph)
        return net


class Assembler:
    def __init__(self, read_path=None):
        self.seqs = None
        self.contigs = None
        self.debrujin_graph = Graph()
        self.K = None
        if read_path is not None:
            self.load_read_data(read_path)

    def load_read_data(self, read_path):
        if not os.path.exists(read_path):
            print(f"File does not exist: {read_path}")
            return -1
        with open(read_path, 'r') as f:
            reads_data = f.readlines()
        self.seqs = []
        for idx, r in enumerate(reads_data):
            if idx % 4 == 1:
                self.seqs.append(r[:-1])

        return 1

    def create_debrujin_graph(self, k):
        self.K = k
        self.debrujin_graph = Graph()
        for s in self.seqs:
            for b in range(len(s)-k):
                v1 = s[b:b+k-1]
                v2 = s[b+1:b+k]
                self.debrujin_graph.add_vertex(v1)
                self.debrujin_graph.add_vertex(v2)
                self.debrujin_graph.add_edge(v1, v2)
        print(f"Create a debrujin graph with {len(self.debrujin_graph.get_vertices())} vertices")

    def plot_debrujin_graph(self):
        net = self.debrujin_graph.plot()
        net.show(f'debrujin_graph_k_{self.K}.html')
    def assemble(self):
        pass

    def calculate_metrics(self):
        pass



