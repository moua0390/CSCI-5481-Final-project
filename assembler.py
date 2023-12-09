import os.path
import os
from pyvis.network import Network
import networkx as nx


class Graph:
    def __init__(self):
        # Initialize an empty dictionary to store the graph
        self.graph = nx.DiGraph()
        self.idx2label = {}

    def add_vertex(self, vertex):
        # Add a vertex to the graph
        if vertex not in self.idx2label.values():
            current_idx = len(self.idx2label.keys())
            self.graph.add_node(current_idx, label=vertex)
            self.idx2label[current_idx] = vertex

    def add_edge(self, vertex1, vertex2):
        # Add an edge between two vertices
        idx_list = list(self.idx2label.keys())
        label_list = list(self.idx2label.values())
        v1_idx = idx_list[label_list.index(vertex1)]
        v2_idx = idx_list[label_list.index(vertex2)]
        self.graph.add_edge(v1_idx, v2_idx)

    def calculate_degree(self, v):
        return {'in': self.graph.in_degree[v], 'out': self.graph.out_degree[v]}

    def get_vertices(self):
        # Get all vertices in the graph
        return self.graph.nodes

    def get_edges(self):
        # Get all edges in the graph
        return self.graph.edges(data=False)

    def remove_nodes(self, v_list):
        for v in v_list:
            print(f'Vertex: {v} which is correspond to k-1 mer {self.idx2label[v]} has been removed in tip removal procedure')
            self.graph.remove_node(v)
            del self.idx2label[v]

    def get_nodes_neighbor(self, v):
        return {'out': [out_v for _, out_v in self.graph.out_edges(v)],
                'in': [in_v for in_v, _ in self.graph.in_edges(v)]}

    def merge_vertices_path(self, v_path, k):
        print(v_path)
        v_source = v_path[0]
        for v_dest in v_path[1:]:
            print(v_source, v_dest)
            nx.contracted_nodes(self.graph, v_source, v_dest, self_loops=False, copy=False)
            del self.graph.nodes[v_source]['contraction']
            v_source_label = self.idx2label[v_source]
            v_dest_label = self.idx2label[v_dest]

            new_label =  v_source_label + v_dest_label[-1]
            attr = {v_source: {'label': new_label}}
            nx.set_node_attributes(self.graph, attr)

            self.idx2label[v_source] = new_label
            del self.idx2label[v_dest]

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

    def plot_debrujin_graph(self, name='debrujin_graph'):
        net = self.debrujin_graph.plot()
        net.show(f'{name}_k_{self.K}.html')

    def simplify_debrujin_graph(self):
        all_edges = self.debrujin_graph.get_edges()
        simple_edges = {}
        simple_paths = []
        reversed_simple_edges = {}
        for v1, v2 in all_edges:
            if self.debrujin_graph.calculate_degree(v1)['out'] == 1 and self.debrujin_graph.calculate_degree(v2)['in'] == 1:
                simple_edges[v1] = v2
                reversed_simple_edges[v2] = v1

        for v in simple_edges.keys():
            if v in [y for x in simple_paths for y in x[:-1]]:
                continue
            simple_paths.append([v])
            v_cur = simple_edges[v]
            while True:
                simple_paths[-1].append(v_cur)
                if v_cur in simple_edges.keys():
                    v_cur = simple_edges[v_cur]
                else:
                    break
            if v in reversed_simple_edges.keys():
                v_cur = reversed_simple_edges[v]
                while True:
                    simple_paths[-1].insert(0, v_cur)
                    if v_cur in reversed_simple_edges.keys():
                        v_cur = reversed_simple_edges[v_cur]
                    else:
                        break

        for s in simple_paths:
            self.debrujin_graph.merge_vertices_path(s, self.K)

    def _remove_tips(self):
        all_vertex = self.debrujin_graph.get_vertices()
        tip_nodes = []
        for v in all_vertex:
            degree = self.debrujin_graph.calculate_degree(v)
            if degree['in'] == 0 or degree['out'] == 0:
                tip_nodes.append(v)
        self.debrujin_graph.remove_nodes(tip_nodes)

    def correct_debrujin_graph_errors(self):
        self._remove_tips()

    def assemble(self):
        all_vertex = list(self.debrujin_graph.get_vertices())
        # Subfunction for determining the best node to start at
        def start_node():
            start = 0
            for v in all_vertex:
                degree = self.debrujin_graph.calculate_degree(v)
                if degree['out'] - degree['in'] == 1:
                    return v
                if degree['out'] > 0:
                    start = v
            return start
        v = start_node()
        eulerian_path = [v]
        neighbor = self.debrujin_graph.get_nodes_neighbor(v)
        while neighbor['out']:
            # Ensure children nodes are in a list to work properly with for-loop
            if not isinstance(neighbor['out'], list):
                neighbor['out'] = [neighbor['out']]
            # Sort children nodes in descending order by the length of their sequence.
            # That way, the longest sequence is chosen to be a part of the Eulerian path.
            neighbor['out'].sort(key=lambda o: len(self.debrujin_graph.idx2label[o]), reverse=True)
            for i, n in enumerate(neighbor['out'], 1):
                if n not in eulerian_path:
                    v = n
                    eulerian_path.append(v)
                    neighbor = self.debrujin_graph.get_nodes_neighbor(v)
                    break
                # If all children nodes are already in path, set outgoing nodes to None to end the loop on the next iteration
                if i == len(neighbor['out']):
                    neighbor['out'] = None
        dropped_nodes = list(set(all_vertex).difference(eulerian_path))
        self.debrujin_graph.remove_nodes(dropped_nodes)
        print(f'The following {len(dropped_nodes)} nodes have been dropped to generate a eulerized graph: {dropped_nodes}')
        assembled_genome = ''
        for e in eulerian_path:
            child_seq = self.debrujin_graph.idx2label[e]
            if assembled_genome == '':
                assembled_genome = child_seq
            else:
                child_seq_slice = child_seq
                child_idx = assembled_genome.find(child_seq)
                while child_idx == -1:
                    child_seq_slice = child_seq_slice[:-1]
                    child_idx = assembled_genome.find(child_seq_slice)
                assembled_genome = assembled_genome[:child_idx] + child_seq
        print(f'Assembled genome: {assembled_genome}')
        return assembled_genome

    def calculate_metrics(self):
        pass




# simple_path = [v1, v2]
# ###
# v_cur = v2
# while True:
#     if self.debrujin_graph.calculate_degree(v_cur)['out'] == 1:
#         v_next = self.debrujin_graph.get_nodes_neighbor(v_cur)['out'][0]
#         if self.debrujin_graph.calculate_degree(v_next)['in'] == 1:
#             simple_path.append(v_next)
#         else:
#             break
#         v_cur = v_next
#     else:
#         break
#
# v_cur = v1
# while True:
#     if self.debrujin_graph.calculate_degree(v_cur)['in'] == 1:
#         v_prev = self.debrujin_graph.get_nodes_neighbor(v_cur)['in'][0]
#         if self.debrujin_graph.calculate_degree(v_prev)['out'] == 1:
#             simple_path.index(0, v_prev)
#         else:
#             break
#         v_cur = v_prev
#     else:
#         break
