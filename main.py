from assembler import Assembler



test_assembler = Assembler()
test_assembler.load_read_data('./lighter_out/sars_spike_protein_raw_reads_trimmed.cor.fq')
# test_assembler.load_read_data('./lighter_out/test.fq')
test_assembler.create_debrujin_graph(k=31)
# test_assembler.plot_debrujin_graph('raw')
test_assembler.simplify_debrujin_graph()
test_assembler.plot_debrujin_graph('simplify')
test_assembler.correct_debrujin_graph_errors()
test_assembler.plot_debrujin_graph('pruned')
test_assembler.assemble()
test_assembler.plot_eulerian_path('eulerized')
print(test_assembler.calculate_metrics())