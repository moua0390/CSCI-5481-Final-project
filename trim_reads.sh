mkdir trimmed_reads

read1=sars_spike_protein_raw_reads_1.fastq
read2=sars_spike_protein_raw_reads_2.fastq

output1_paired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_paired.fastq
output1_unpaired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_unpaired.fastq

output2_paired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_paired.fastq
output2_unpaired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_unpaired.fastq

logfile=trimmed_reads/log.txt

java -jar trimmomatic/trimmomatic-0.39.jar PE -trimlog $logfile $read1 $read2 \
$output1_paired $output1_unpaired \
$output2_paired $output2_unpaired \
ILLUMINACLIP:trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True \
LEADING:3 TRAILING:3 MINLEN:36
