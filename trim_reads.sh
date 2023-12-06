mkdir -p trimmed_reads

read1=sars_spike_protein_raw_reads_1.fastq
read2=sars_spike_protein_raw_reads_2.fastq

head -n 10000 sars_spike_protein_raw_reads.fastq > $read1
tail -n 10000 sars_spike_protein_raw_reads.fastq > $read2

output1_paired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_paired.fastq
output1_unpaired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_unpaired.fastq

output2_paired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_paired.fastq
output2_unpaired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_unpaired.fastq

logfile=trimmed_reads/log.txt
statsfile=trimmed_reads/stats.txt

java -jar packages/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 \
-trimlog $logfile -summary $statsfile \
$read1 $read2 \
$output1_paired $output1_unpaired \
$output2_paired $output2_unpaired \
SLIDINGWINDOW:4:25 LEADING:35 TRAILING:35 MINLEN:35
