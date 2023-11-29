mkdir -p fastqc_output

read1=sars_spike_protein_raw_reads_1.fastq.gz
read2=sars_spike_protein_raw_reads_2.fastq.gz

trimmed_read1=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_paired.fastq.gz
trimmed_read2=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_paired.fastq.gz

packages/FastQC/fastqc $read1 $read2 $trimmed_read1 $trimmed_read2 -o fastqc_output
