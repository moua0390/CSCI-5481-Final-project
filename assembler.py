import os.path
import os
from pyvis.network import Network
import networkx as nx


class Graph:
    def __init__(self):
        # Initialize an empty dictionary to store the graph
        # self.graph = nx.MultiDiGraph()
        self.graph = nx.DiGraph()
        self.idx2label = {}
        self.label2idx = {}

    def add_vertex(self, vertex):
        # Add a vertex to the graph
        if vertex not in self.idx2label.values():
            current_idx = len(self.idx2label.keys())
            self.graph.add_node(current_idx, label=vertex)
            self.idx2label[current_idx] = vertex
            self.label2idx[vertex] = current_idx

    def add_edge(self, vertex1, vertex2):
        # Add an edge between two vertices
        v1_idx = self.label2idx[vertex1]
        v2_idx = self.label2idx[vertex2]
        if self.graph.has_edge(v1_idx, v2_idx):
            cur_weight = self.graph[v1_idx][v2_idx]['w']
            self.graph[v1_idx][v2_idx]['w'] = cur_weight + 1
            self.graph[v1_idx][v2_idx]['label'] = str(cur_weight + 1)
        else:
            self.graph.add_edge(v1_idx, v2_idx, w=1, label='1')

    def calculate_degree(self, v):
        return {'in': self.graph.in_degree(v, weight='w'),
                'out': self.graph.out_degree(v, weight='w'),
                'in_nodes': self.graph.in_degree(v),
                'out_nodes': self.graph.out_degree(v)}

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
            del self.label2idx[self.idx2label[v]]
            del self.idx2label[v]

    def get_nodes_neighbor(self, v):
        return {'out': [out_v for _, out_v in self.graph.out_edges(v)],
                'in': [in_v for in_v, _ in self.graph.in_edges(v)]}

    def merge_vertices_path(self, v_path):
        print("Path is merging: ", v_path)
        v_source = v_path[0]
        path_weight_list = []
        path_weight = self.graph.out_degree(v_source, weight='w')
        for v_dest in v_path[1:]:
            if self.graph.in_degree(v_dest, weight='w') < path_weight:
                path_weight = self.graph.in_degree(v_dest, weight='w')

            path_weight_list.append(self.graph.in_degree(v_dest, weight='w'))
            nx.contracted_nodes(self.graph, v_source, v_dest, self_loops=False, copy=False)
            del self.graph.nodes[v_source]['contraction']
            v_source_label = self.idx2label[v_source]
            v_dest_label = self.idx2label[v_dest]

            new_label = v_source_label + v_dest_label[-1]
            attr = {v_source: {'label': new_label}}
            nx.set_node_attributes(self.graph, attr)

            self.idx2label[v_source] = new_label
            self.label2idx[new_label] = v_source
            del self.idx2label[v_dest]
            del self.label2idx[v_source_label]
            del self.label2idx[v_dest_label]

        print("Path Weight: ", path_weight_list)
        path_end_branches = list(self.graph.out_edges(v_source))
        if len(path_end_branches) > 0:
            next_v = path_end_branches[0][1]
            self.graph[v_source][next_v]['w'] = path_weight
            self.graph[v_source][next_v]['label'] = str(path_weight)

    def plot(self, width=1000, height=1000):
        copy_graph = self.graph.copy()
        net = Network(f'{width}px', f'{height}px', directed=True)
        net.from_nx(copy_graph)
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
        seq_lengths = []
        for idx, r in enumerate(reads_data):
            if idx % 4 == 1:
                self.seqs.append(r[:-1])
                seq_lengths.append(len(r[:-1]))

        print("Loading Data..")
        print(f"Min Length: {min(seq_lengths)}, Max Length: {max(seq_lengths)},"
              f" Average: {sum(seq_lengths)/len(seq_lengths)}")

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
        print(f"Create a debrujin graph with {len(self.debrujin_graph.get_vertices())} vertices and"
              f" {len(self.debrujin_graph.get_edges())} edges")

    def plot_debrujin_graph(self, name='debrujin_graph'):
        net = self.debrujin_graph.plot()
        net.show(f'{name}_k_{self.K}.html')

    def simplify_debrujin_graph(self):
        all_edges = self.debrujin_graph.get_edges()
        simple_edges = {}
        simple_paths = []
        reversed_simple_edges = {}
        for v1, v2 in all_edges:
            if (self.debrujin_graph.calculate_degree(v1)['out_nodes'] == 1 and
                    self.debrujin_graph.calculate_degree(v2)['in_nodes'] == 1):
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
            if 4073 in s:
                print(1)
            self.debrujin_graph.merge_vertices_path(s)

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
        pass

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