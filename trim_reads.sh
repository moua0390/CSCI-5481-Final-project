mkdir -p trimmed_reads

input=sars_spike_protein_raw_reads.fastq
output=trimmed_reads/sars_spike_protein_raw_reads_trimmed.fastq

logfile=trimmed_reads/log.txt
statsfile=trimmed_reads/stats.txt

java -jar packages/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -phred33 \
-trimlog $logfile -summary $statsfile $input $output \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:35

mkdir -p lighter_out
Lighter_exec/lighter -r trimmed_reads/sars_spike_protein_raw_reads_trimmed.fastq -K 5 5000 -od lighter_out/
