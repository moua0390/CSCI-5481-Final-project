from assembler import Assembler
from Bio import Align
import matplotlib.pyplot as plt
from numpy import random
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

# Retrieve the first optimal alignment between the longest gene and contig
separate_data = read_file('./SARS-CoV-2_separate_genes.fna')
longest_gene_seq = max(separate_data[1], key=len)
longest_gene_id = separate_data[0][separate_data[1].index(longest_gene_seq)]

assembled_data = read_file('./contigs.fna')
longest_contig_seq = max(assembled_data[1], key=len)
longest_contig_id = assembled_data[0][assembled_data[1].index(longest_contig_seq)]

optimal_alignment = None
aligner = Align.PairwiseAligner()
with open('alignments.txt', 'w') as alignment_file:
    alignments = aligner.align(longest_gene_seq, longest_contig_seq)
    alignment_file.write(f'{alignments}')
    optimal_alignment = alignments[0]
    alignment_file.write(f'{longest_gene_id}, {longest_contig_id}\nScore = {optimal_alignment.score}\n{optimal_alignment}\n')

# Generate scatterplot
num_points = len(longest_contig_seq)
x = []
y = [*range(num_points)]
position = 0
with open('alignments.txt', 'r') as alignment_file:
    all_lines = alignment_file.readlines()
    for idx in range(3, len(all_lines), 4):
        for c in all_lines[idx]:
            if c != '|' and c != '-':
                continue
            else:
                position += 1
                if c == '|':
                    x.append(position)
                    if len(x) == num_points:
                        break
colors = random.rand(num_points)
plt.scatter(x, y, c=colors)
plt.title("Base position in longest assembled contig vs. reference gene")
plt.xlabel("Position in longest gene")
plt.ylabel("Position in longest contig")
plt.show()
