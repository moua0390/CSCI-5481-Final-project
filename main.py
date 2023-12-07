from assembler import Assembler


test_assembler = Assembler()
test_assembler.load_read_data('./lighter_out/sars_spike_protein_raw_reads_1_trimmed_paired.cor.fq')
test_assembler.create_debrujin_graph(k=23)
test_assembler.plot_debrujin_graph()
