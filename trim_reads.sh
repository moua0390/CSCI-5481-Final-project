mkdir trimmed_reads

read1=sars_spike_protein_raw_reads_1.fq
read2=sars_spike_protein_raw_reads_2.fq

head -n 10000 sars_spike_protein_raw_reads.fastq > $read1
tail -n 10000 sars_spike_protein_raw_reads.fastq > $read2

gzip $read1
gzip $read2

output1_paired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_paired.fq.gz
output1_unpaired=trimmed_reads/sars_spike_protein_raw_reads_1_trimmed_unpaired.fq.gz

output2_paired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_paired.fq.gz
output2_unpaired=trimmed_reads/sars_spike_protein_raw_reads_2_trimmed_unpaired.fq.gz

logfile=trimmed_reads/log.txt

java -jar trimmomatic/trimmomatic-0.39.jar PE -trimlog $logfile $read1.gz $read2.gz \
$output1_paired $output1_unpaired \
$output2_paired $output2_unpaired \
ILLUMINACLIP:trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:True \
LEADING:3 TRAILING:3 MINLEN:36
