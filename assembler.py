import os.path
import os
from pyvis.network import Network
import networkx as nx
import numpy as np


class Graph:
    def __init__(self, graph=None, idx2label_dict=None, label2idx_dict=None):
        # Initialize an empty dictionary to store the graph
        self.graph = nx.DiGraph() if graph is None else graph
        self.idx2label = {} if idx2label_dict is None else idx2label_dict
        self.label2idx = {} if label2idx_dict is None else label2idx_dict

    def __copy__(self):
        return Graph(self.graph.copy(), self.idx2label.copy(), self.label2idx.copy())

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
        self.contigs = []
        self.assembled_contigs = []
        self.debrujin_graph = Graph()
        self.eulerian_paths = Graph()
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

    def plot_eulerian_path(self, name='eulerian_graph'):
        net = self.eulerian_paths.plot(width=400, height=400)
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
        all_vertex = list(self.debrujin_graph.get_vertices())
        # Subfunction for determining the best node to start at
        def start_node():
            start = all_vertex[0]
            for v in all_vertex:
                degree = self.debrujin_graph.calculate_degree(v)
                if degree['out'] - degree['in'] == 1:
                    return v
                if degree['out'] > 0:
                    start = v
            return start
        start_nodes = []
        visited = []
        # try:
        while len(all_vertex) > 0:
            v = start_node()
            start_nodes.append(v)
            eulerian_path = [v]
            neighbor = self.debrujin_graph.get_nodes_neighbor(v)
            while neighbor['out']:
                # Ensure children nodes are in a list to work properly with for-loop
                if not isinstance(neighbor['out'], list):
                    neighbor['out'] = [neighbor['out']]
                # Sort children nodes in descending order by the length of their sequence.
                # That way, the longest sequence is chosen to be a part of the Eulerian path.
                neighbor['out'].sort(key=lambda node: len(self.debrujin_graph.idx2label[node]), reverse=True)
                count = 0
                for n in neighbor['out']:
                    if n not in visited:
                        v = n
                        eulerian_path.append(v)
                        visited.append(v)
                        neighbor = self.debrujin_graph.get_nodes_neighbor(v)
                        break
                    count += 1
                # Check if all children nodes are already in path
                if count == len(neighbor['out']):
                    break
            print("PATH: ", eulerian_path)
            self.contigs.append(eulerian_path)
            all_vertex = list(set(all_vertex).difference(eulerian_path))
        print("Assembly Done........................")

        sorted(self.contigs, key=lambda x: len(x), reverse=True)
        self.contigs = [c for c in self.contigs if len(c) > 1]

        self.eulerian_paths = self.debrujin_graph.__copy__()
        for idx, contig in enumerate(self.contigs, 1):
            contig_color = (int(np.random.random()*255),  int(np.random.random()*255), int(np.random.random()*255))
            assembled = ''
            for ii, c in enumerate(contig):
                if ii != len(contig)-1:
                    next_node = contig[ii+1]
                    self.eulerian_paths.graph[c][next_node]['color'] =\
                        f'#{hex(contig_color[0])[2:]}{hex(contig_color[1])[2:]}{hex(contig_color[2])[2:]}'
                    self.eulerian_paths.graph[c][next_node]['value'] = 15

                child_seq = self.eulerian_paths.idx2label[c]
                if assembled == '':
                    assembled = child_seq
                else:
                    assembled = assembled + child_seq[self.K-2:]
            self.assembled_contigs.append(assembled)
            print(f'Assembled contig {idx}: {assembled}')

    def calculate_metrics(self):
        total_length = sum([len(c) for c in self.assembled_contigs])
        N50, L50, N90 = None, None, None
        cur_sum = 0
        for idx, c in enumerate(self.assembled_contigs, 1):
            cur_sum += len(c)
            if cur_sum >= 0.5 * total_length and N50 is None:
                N50 = len(c)
                L50 = idx

            if cur_sum >= 0.9 * total_length and N90 is None:
                N90 = len(c)

        return {'N50': N50, 'L50': L50, 'N90': N90}

