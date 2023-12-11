from assembler import Assembler
from Bio import Align
import matplotlib.pyplot as plt
from numpy import random
import os.path

def read_contigs_file(file_path):
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
            ids.append(r[1:-1])
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

# Retrieve sequence from whole genome file
whole_genome_id = 'SARS-CoV-2_whole_genome'
whole_genome_seq = None
with open(f'{whole_genome_id}.fna', 'r') as whole_genome_file:
    whole_genome_seq = whole_genome_file.readlines()[1]

# Retrieve longest contig
assembled_data = read_contigs_file('./contigs.fna')
all_contigs_id = assembled_data[0]
all_contigs_seq = assembled_data[1]
longest_contig_seq = max(assembled_data[1], key=len)
longest_contig_id = assembled_data[0][assembled_data[1].index(longest_contig_seq)]

optimal_alignment = None
aligner = Align.PairwiseAligner()
aligner.match_score = 1
aligner.mismatch_score = -2
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
with open('alignments.txt', 'w') as alignment_file:
    # Align all contigs against whole genome sequence
    for i in range(len(all_contigs_seq)):
        contig_id = all_contigs_id[i]
        contig_seq = all_contigs_seq[i]
        alignments = aligner.align(whole_genome_seq, contig_seq)
        # Retrieve only the first optimal alignment to save computation time
        first_alignment = alignments[0]
        alignment_file.write(f'{whole_genome_id}, {contig_id}\nScore = {first_alignment.score}\n{first_alignment}\n')

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
plt.scatter(x, y, s=1, c=colors)
plt.title("Base position in longest assembled contig vs. reference gene")
plt.xlabel("Position in longest gene")
plt.ylabel("Position in longest contig")
plt.show()
