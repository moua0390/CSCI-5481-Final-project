from assembler import Assembler
from pyvis.network import Network
#
# import networkx as nx
# net = nx.DiGraph()
# net.add_node(1, label='s')
# net.add_node(2, label='p')
# net.add_node(3, label='t')
# net.add_node(4, label='h')
# net.add_node(5, label='aa')
# net.add_node(6, label='bb')
# net.add_node(7, label='cc')
# net.add_node(8, label='dd')
# net.add_edge(1, 2)
# net.add_edge(2, 3)
# net.add_edge(3, 4)
# net.add_edge(5, 1)
# net.add_edge(6, 1)
# net.add_edge(4, 7)
# net.add_edge(4, 8)


# net2 = nx.contracted_nodes(net, 2, 4, self_loops=False)
# net3 = nx.contracted_nodes(net2, 1, 2, self_loops=False)


# for n in net.edges(data=True):
#     print(n)
#
# nx.contracted_nodes(net, 1, 2, self_loops=False, copy=False)
# print('####')
#
# for n in net.edges():
#     print(n)

# nt = Network('1000px', '1000px', directed=True)
# nt.from_nx(net)
# nt.show('test1.html')

# print(net.out_edges(2))

# for n in net3.nodes(data=True):
#     print(n)
#
# atts = {1:{'label':'adfasdf'}}
# nx.set_node_attributes(net3, atts )
# del net3.nodes[1]['contraction']
#
# print("################")
# for n in net3.nodes(data=True):
#     print(n)
#
# for v1, v2, d in net3.edges(data=True):
#     if 'contraction' in d.keys():
#         del net3[v1][v2]['contraction']


test_assembler = Assembler()
test_assembler.load_read_data('./lighter_out/sars_spike_protein_raw_reads_trimmed.cor.fq')
# test_assembler.load_read_data('./lighter_out/test.fq')
test_assembler.create_debrujin_graph(k=31)
# test_assembler.plot_debrujin_graph('raw')
test_assembler.simplify_debrujin_graph()
test_assembler.plot_debrujin_graph('simplify')
# test_assembler.correct_debrujin_graph_errors()
# test_assembler.plot_debrujin_graph('pruned')
test_assembler.assemble()
# test_assembler.plot_debrujin_graph('eulerized')
