from assembler import Assembler
from Bio import Align
import matplotlib.pyplot as plt
from numpy import random
import os.path

def alignment():
    # Retrieve sequence from whole genome file
    whole_genome_id = 'sars_spike_protein_assembled'
    whole_genome_seq = None
    with open(f'{whole_genome_id}.fna', 'r') as whole_genome_file:
        whole_genome_seq = whole_genome_file.readlines()[1]
    
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.mode = 'local'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    with open('alignments.txt', 'w') as alignment_file:
        # Align all contigs against whole genome sequence
        for key in (contig_dict):
            contig_id = key
            contig_seq = contig_dict[key]
            alignments = aligner.align(whole_genome_seq, contig_seq)
            # Retrieve only the first optimal alignment to save computation time
            first_alignment = alignments[0]
            alignment_file.write(f'{whole_genome_id} <- {contig_id}, Score = {first_alignment.score}\n{first_alignment}\n')

def create_scatterplot(contig_id):
    num_bases = len(contig_dict[contig_id])
    x = []
    y = [*range(num_bases)]
    position = 0
    with open('alignments.txt', 'r') as alignment_file:
        for line in alignment_file:
            if contig_id in line:
                break
            position += 1
        for curr_line in alignment_file:
            if 'target' in curr_line or 'query' in curr_line:
                continue
            if 'Contig' in curr_line:
                break
            for c in curr_line:
                if c != '|' and c != '-' and '.':
                    continue
                else:
                    if c == '|':
                        x.append(position)
                        if len(x) == num_bases:
                            break
                        position += 1
            if len(x) == num_bases:
                break
    while len(x) > len(y):
        x.pop()
    while len(x) < len(y):
        y.pop()
    colors = random.rand(len(x))
    plt.scatter(x, y, s=1, c=colors)
    plt.title(f'Base position in {contig_id} vs. reference genome')
    plt.xlabel('Position in ref genome')
    plt.ylabel('Position in contig')
    plt.show()

def read_contigs_file(file_path):
    if not os.path.exists(file_path):
        print(f'File does not exist: {file_path}')
    with open(file_path, 'r') as f:
        reads_data = f.readlines()
    sequences = []
    ids = []
    for idx, r in enumerate(reads_data):
        if idx % 2 != 0:
            sequences.append(r[:-1])
        else:
            ids.append(r[1:-1])
    return dict(zip(ids, sequences))



test_assembler = Assembler()
test_assembler.load_read_data('./lighter_out/sars_spike_protein_raw_reads_trimmed.cor.fq')
test_assembler.create_debrujin_graph(k=31)
test_assembler.plot_debrujin_graph('raw')
test_assembler.simplify_debrujin_graph()
test_assembler.plot_debrujin_graph('simplify')
test_assembler.correct_debrujin_graph_errors()
test_assembler.plot_debrujin_graph('pruned')
test_assembler.assemble()
test_assembler.plot_eulerian_path('eulerized')
print(test_assembler.calculate_metrics())

# Save all contigs in a dictionary
contig_dict = read_contigs_file('./contigs.fna')
# Align all contigs to whole genome sequence. Save alignments to a text file
alignment()
# Generate a scatterplot of the base positions in the whole genome file and given contig
create_scatterplot('Contig 1')
create_scatterplot('Contig 2')
create_scatterplot('Contig 3')
