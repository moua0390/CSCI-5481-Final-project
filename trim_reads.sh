mkdir -p trimmed_reads

input=sars_spike_protein_raw_reads.fastq
output=trimmed_reads/sars_spike_protein_raw_reads_trimmed.fastq

logfile=trimmed_reads/log.txt
statsfile=trimmed_reads/stats.txt

java -jar packages/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -phred33 \
-trimlog $logfile -summary $statsfile $input $output \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:35

mkdir -p lighter_out
packages/Lighter/lighter -r trimmed_reads/sars_spike_protein_raw_reads_trimmed.fastq -k 17 5000 0.035 -t 10 -od lighter_out/
