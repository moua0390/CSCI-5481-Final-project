mkdir -p fastqc_output

raw_reads=sars_spike_protein_raw_reads.fastq
trimmed_reads=trimmed_reads/sars_spike_protein_raw_reads_trimmed.fastq
corrected_reads=lighter_out/sars_spike_protein_raw_reads_trimmed.cor.fq

packages/FastQC/fastqc $raw_reads $trimmed_reads $corrected_reads -o fastqc_output
