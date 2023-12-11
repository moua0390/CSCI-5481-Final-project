from assembler import Assembler
from Bio import Align, AlignIO
import os.path

def read_file(file_path):
    if not os.path.exists(file_path):
        print(f"File does not exist: {file_path}")
    with open(file_path, 'r') as f:
        reads_data = f.readlines()
    sequences = []
    ids = []
    for idx, r in enumerate(reads_data):
        if idx % 2 != 0:
            sequences.append(r[:-1])
        else:
            if '>lcl|' in r:
                id = r[:r.index(' ')].strip('>lcl|')
            else:
                id = r[1:-1]
            ids.append(id)
    return (ids, sequences)

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
test_assembler.plot_eulerian_path('eulerized')
print(test_assembler.calculate_metrics())

separate_data = read_file('./SARS-CoV-2_separate_genes.fna')
assembled_data = read_file('./contigs.fna')
longest_contig_seq = max(assembled_data[1], key=len)
longest_contig_id = assembled_data[0][assembled_data[1].index(longest_contig_seq)]
aligner = Align.PairwiseAligner()
with open('alignments.txt', 'w') as alignment_file:
    for idx, s in enumerate(separate_data[1]):
        alignments = aligner.align(s, longest_contig_seq)
        alignment_file.write(f'{alignments}')
        for alignment in alignments:
            alignment_file.write(f'{separate_data[0][idx]}, {longest_contig_id}\nScore = {alignment.score}\n{alignment}\n')
