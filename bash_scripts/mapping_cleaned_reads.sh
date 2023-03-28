#!/bin/bash

{
usage="$(basename "$0") [-h] [-I <SRA_list_Illumina>] [-N <SRA_list_Nanopore>] [-g <reference_genome>] [-d <working_directory]
This program will trim and map reads from given SRA sequences files to compare mapping of different NGS datasets.
    -h  show this help text
    -I  File or path to SRA accession list for Illumina data in tabular format
    -N  File or path to SRA accession list for Nanopore data in tabular format
    -g  The SARS-CoV-2 reference genome (in fasta format). It is assumed that the genome is within the References directory in -d 
    -d  Path to the working directory (the main directory for the repository)
    -t  Number of CPU processors"
options=':h:N:I:g:t:d:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    I) I=$OPTARG;;
    N) N=$OPTARG;;
    g) g=$OPTARG;;
    t) t=$OPTARG;;
    d) d=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

echo ""
echo "Working directory: " $d
echo "Illumina data:     " $I
echo "Nanopore data:     " $N
echo "Genome file:       " $g
echo "Threads:           " $t
echo ""

# mandatory arguments
if [ ! "$I" ] || [ ! "$N" ] || [ ! "$g" ] || [ ! "$t" ] || [ ! "$d" ]; then
  echo "arguments -I, -N, -g, -t, and -d  must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`


# Setup working environment
## -- Illumina
mkdir -p Illumina/sra_files
mkdir -p Illumina/raw_reads
mkdir -p Illumina/cleaned_reads
mkdir -p Illumina/mapped_reads
## -- Nanopore
mkdir -p Nanopore/sra_files
mkdir -p Nanopore/raw_reads
mkdir -p Nanopore/cleaned_reads
mkdir -p Nanopore/mapped_reads

##################################################################################
#  Download fastq files for 22 samples from the Heikema et. al (2020) paper      # 
##################################################################################

module load sra-toolkit/3.0.2

cd ${d}
echo "Downloading Illumina SRA files from the given list of accessions"
cd Illumina/sra_files
prefetch --max-size 800G -O ./ --option-file ${d}/${I}
echo "Converting Illumina SRA files to fastq.gz"
ls -p | grep ERR > sra_dirs
while read i; do mv "$i"*.sra .; rmdir "$i"; done<sra_dirs
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done
mv *.fastq.gz ../raw_reads

cd ${d}
pwd
echo "Downloading Nanopore SRA files from the given list of accessions"
cd Nanopore/sra_files
prefetch --max-size 800G -O ./ --option-file ${d}/${N}
echo "Converting Nanopore SRA files to fastq.gz"
ls -p | grep ERR > sra_dirs
while read i; do mv "$i"*.sra .; rmdir "$i"; done<sra_dirs
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done
mv *.fastq.gz ../raw_reads


module unload sra-toolkit/3.0.2

echo "SRA files were downloaded and .fastq files extracted in Illumina and Nanopore directories"
echo ""
echo "Done"
echo ""

########################################################################################
# Trimming downloaded Illumina datasets with fastp, using 16 threads (-w option) #
########################################################################################

module load fastp/0.20.1

echo "Trimming downloaded datasets with fastp."
echo ""

function trim_reads {
cd ${d}/$1/raw_reads
ls *.fastq.gz | cut -d "." -f "1" | sort > fastq_list
while read sample; do fastp -w ${t} -i ${sample}.fastq.gz -o ../cleaned_reads/${sample}_cleaned.fastq.gz
done<fastq_list
cd ..
}

trim_reads Illumina
trim_reads Nanopore

  
module unload fastp/0.20.1
    
######################
# Indexing Reference #
######################

module load bwa/2020_03_19

echo "Indexing Reference"
echo ""
cd ${d}
bwa index ${d}/Reference/${g}

module unload bwa/2020_03_19

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

module load bwa/2020_03_19
module load samtools/1.16

echo "Aligning Illumina datasets againts reference with bwa mem, using $t threads."
echo ""

function align_reads {
cd $1/cleaned_reads
ls *_cleaned.fastq.gz | awk -F'[._]' '{print $1}' | sort > cleaned_list
while read sample; do bwa mem ${d}/Reference/${g} ${sample}_cleaned.fastq.gz > ${d}/$1/mapped_reads/${sample}_mapped.sam
samtools sort ${d}/$1/mapped_reads/${sample}_mapped.sam > ${d}/$1/mapped_reads/${sample}_sorted.bam -@ ${t}
done<cleaned_list
cd ..
}

align_reads Illumina
align_reads Nanopore

module unload bwa/2020_03_19
module unload samtools/1.16

########## Done ###########

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

mdate=`date +'%d/%m/%Y %H:%M:%S'`
# NOTE: If you are running a mac and having trouble with the code below,
# ----- try using 'vm_stat' instead of 'vmstat'
mcpu=$[100-$(vmstat 1 2|tail -1|awk '{print $15}')]%
mmem=`free | grep Mem | awk '{print $3/$2 * 100.0}'`
echo "$mdate | $mcpu | $mmem" >> ./stats-cpu
###############################################################
#
} | tee logfile
