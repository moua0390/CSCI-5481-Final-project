# CSCI 5481 Final project (Group 22)

## Group members
+ Hmong Moua (`moua0390`)
+ Parsa Moradi (`moradi`)

## Proposal
We will attempt to perform Covid genome assembly (project option 2). This will involve applying quality trimming, correcting sequencing errors, dividing into k-mers (testing various values for k), and constructing as well as plotting De Bruijn graph on raw shotgun sequencing data of the SARS-CoV2 spike protein. Afterwards, we will attempt to assemble the genome by finding an Euclearian path, handling arising issues then align an assembled contig to the spike protein from the separate-genes file in Homework 1 and create a scatter plot to compare the location of their bases. Finally we will compare our assembly performance via known metrics (such as N50) with benchmark methods.

## Approach
### Step 1: Perform quality trimming/filtering using Trimmomatic
+ Quality trimming and filtering is performed on raw shotgun sequencing data for the SARS-CoV2 spike protein (`sars_spike_protein_raw_reads.fastq`) using Trimmomatic. This can be performed by running the commands in `trim_reads.sh`. To provide some clarification, here is a breakdown of some of the commands in the file:
  + We identified that the raw sequencing data is paired-end, so the following commands divide the first 10,000 lines and last 10,000 lines, which represent the forward and reverse reads, into two separate files.
    ```
    read1=sars_spike_protein_raw_reads_1.fastq
    read2=sars_spike_protein_raw_reads_2.fastq
    
    head -n 10000 sars_spike_protein_raw_reads.fastq > $read1
    tail -n 10000 sars_spike_protein_raw_reads.fastq > $read2
    ```
  + Here is the actual execution of Trimmomatic.
    ```
    java -jar packages/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 \
    -trimlog $logfile -summary $statsfile \
    $read1 $read2 \
    $output1_paired $output1_unpaired \
    $output2_paired $output2_unpaired \
    SLIDINGWINDOW:4:25 LEADING:35 TRAILING:35 MINLEN:35
    ```
   + We chose the parameters and their values by analyzing the FastQC output of the data before and after trimming. This can be performed by running the commands in `run_fastqc.sh`. To be specific, the files that were performed on are:
      + `sars_spike_protein_raw_reads_1.fastq`
      + `sars_spike_protein_raw_reads_2.fastq`
      + `sars_spike_protein_raw_reads_1_trimmed_paired.fastq`
      + `sars_spike_protein_raw_reads_2_trimmed_paired.fastq`

    #### Findings before trimming
    + *Adapter Content* showed that the reads did not contain any adapters, so we concluded that they may have been removed beforehand. Thus, we left out the `ILLUMINACLIP` parameter for removing adapters.

      <img src=images/adapter_content.png width=50% height=50%>
      
    + *Per base sequence quality* showed the quality of the bases were lower at the ends. This pushed us to use the `LEADING` and `TRAILING` parameters to improve sequence quality.
      + sars_spike_protein_raw_reads_1.fastq
      <img src=images/reads1_untrimmed_per_base.png width=50% height=50%>
 
      + sars_spike_protein_raw_reads_2.fastq
      <img src=images/reads2_untrimmed_per_base.png width=50% height=50%>
      
    + *Per sequence quality scores* showed that the most common average quality score was 36. This pushed us to use the `MINLEN` parameter and 35 as the threshold for this and aforementioned parameters.
      + sars_spike_protein_raw_reads_1.fastq
      <img src=images/reads1_untrimmed_per_seq.png width=50% height=50%>
 
      + sars_spike_protein_raw_reads_2.fastq
      <img src=images/reads2_untrimmed_per_seq.png width=50% height=50%>
      
    + For further quality filtering, we used the `SLIDINGWINDOW` parameter with a window size of 4 and a threshold of 25 to be more lenient to avoid discarding too much data.

    #### Findings after trimming
    + *Per base sequence quality* showed the quality of the bases improved at the ends after trimming.
      + sars_spike_protein_raw_reads_1_trimmed_paired.fastq
      <img src=images/reads1_trimmed_per_base.png width=50% height=50%>

      + sars_spike_protein_raw_reads_2_trimmed_paired.fastq
      <img src=images/reads2_trimmed_per_base.png width=50% height=50%>
      
    + *Per sequence quality scores* showed that the most common average quality score after trimming was 36 for read 1 and 35 for read 2. This change for read 2 likely due to sequences being dropped.
      + sars_spike_protein_raw_reads_1_trimmed_paired.fastq
      <img src=images/reads1_trimmed_per_seq.png width=50% height=50%>

      + sars_spike_protein_raw_reads_2_trimmed_paired.fastq
      <img src=images/reads2_trimmed_per_seq.png width=50% height=50%>
