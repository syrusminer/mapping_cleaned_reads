# Mapping Cleaned Reads Worksheet

<!--- Write name below --->
## Name: Syrus Miner 

<!--- For this worksheet, answer the following questions --->

## Q1: What does it mean to map/align reads to a reference?

## Answer:
Mapping and alligning reads to a reference genome involves determining the most likely origin of each read to 
the reference sequence. By doing this one can determine genetic variations in the genetic code.   

## Q2: What read mapper does the mapping_cleaned_reads.sh script use?
## Answer: 
bwa

## Q3: Both Illumina and Nanopore reads are used for this assignment. What are the major differences between the methodology used for these sequencing platforms?
## Answer: 
Illumina uses short reads which range from 50 to 300 bp in length while Nanopore uses generates long reads which
can be thousands of basepairs in length. The way in which the data collected for use also varies between the platforms.

## Q4: What differences do you notice between the Illumina and Nanopore raw_data fastq file sizes? Which are larger?
## Answer: 
The Nanopore raw_data files are larger than the Illumina raw_data files. This is due to the fact that
Nanopore reads are longer than the Illumina reads.

## Q5: What differences do you notice between the Illumina and Nanopore cleaned_reads fastq file sizes? Which are larger?
## Answer: 
The opposite is seen compared to question four. Nanopore file sizes for cleaned reads are smaller than the file
sizes for Illumina. This could be due to the fact that Illumina sequencing typically has less errors than Nanopore sequencing.
This is not true always, there seems to be much more variety in the Nanopore files than that of Illumina.

## Q6: What explains the difference in your responses of Q4 and Q5? (HINT: Take a glimpse at the raw data .fastq files themselves)
## Answer: 
This is what I was explaining in question 5, there are more errors made in the Nanopore sequencing than compared to the
Ilumina data, because of this there is a lot more variety in the file sizes of Nanopore data.

## Q7: What is the average read depth for the Illumina data across all samples for the genomic regions that were mapped to?
## Answer: 
37901.4

## Q8: What is the average read depth for the Illumina data across all samples for all genomic regions?
## Answer: 
337.961

## Q9: What is the average read depth for the Nanopore data across all samples for the genomic regions that were mapped to?
## Answer: 
0.0465426

## Q10: What is the average read depth for the Nanopore data across all samples for all genomic regions?
## Answer: 
0.0176322
