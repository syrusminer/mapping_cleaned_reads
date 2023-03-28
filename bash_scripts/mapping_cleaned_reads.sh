#!/bin/bash

{
usage="$(basename "$0") [-h] [-I <SRA_list_Illumina>] [-N <SRA_list_Nanopore>] [-g <reference_genome>] [-d <working_directory]
This program will trim and map reads from given SRA sequences files to compare mapping of different NGS datasets.
    -h  show this help text
    -I  File or path to SRA accession list for Illumina data in tabular format
    -N  File or path to SRA accession list for Nanopore data in tabular format
    -g  The SARS-CoV-2 reference genome (in fasta format). If the genome is located in the working folder, just specify the name.
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

#+ COMPLETED   module load sra-toolkit/3.0.2
#+ COMPLETED   
#+ COMPLETED   cd ${d}
#+ COMPLETED   echo "Downloading Illumina SRA files from the given list of accessions"
#+ COMPLETED   cd Illumina/sra_files
#+ COMPLETED   prefetch --max-size 800G -O ./ --option-file ${d}/${I}
#+ COMPLETED   echo "Converting Illumina SRA files to fastq.gz"
#+ COMPLETED   ls -p | grep ERR > sra_dirs
#+ COMPLETED   while read i; do mv "$i"*.sra .; rmdir "$i"; done<sra_dirs
#+ COMPLETED   SRA= ls -1 *.sra
#+ COMPLETED   for SRA in *.sra; do fastq-dump --gzip ${SRA}
#+ COMPLETED   done
#+ COMPLETED   mv *.fastq.gz ../raw_reads
#+ COMPLETED   
#+ COMPLETED   cd ${d}
#+ COMPLETED   pwd
#+ COMPLETED   echo "Downloading Nanopore SRA files from the given list of accessions"
#+ COMPLETED   cd Nanopore/sra_files
#+ COMPLETED   prefetch --max-size 800G -O ./ --option-file ${d}/${N}
#+ COMPLETED   echo "Converting Nanopore SRA files to fastq.gz"
#+ COMPLETED   ls -p | grep ERR > sra_dirs
#+ COMPLETED   while read i; do mv "$i"*.sra .; rmdir "$i"; done<sra_dirs
#+ COMPLETED   SRA= ls -1 *.sra
#+ COMPLETED   for SRA in *.sra; do fastq-dump --gzip ${SRA}
#+ COMPLETED   done
#+ COMPLETED   mv *.fastq.gz ../raw_reads
#+ COMPLETED   
#+ COMPLETED    
#+ COMPLETED    module unload sra-toolkit/3.0.2
#+ COMPLETED    
#+ COMPLETED    echo "SRA files were downloaded and .fastq files extracted in Illumina and Nanopore directories"
#+ COMPLETED    echo ""
#+ COMPLETED    echo "Done"
#+ COMPLETED    echo ""
#+ COMPLETED    
#+ COMPLETED    ########################################################################################
#+ COMPLETED    # Trimming downloaded Illumina datasets with fastp, using 16 threads (-w option) #
#+ COMPLETED    ########################################################################################
#+ COMPLETED    
#+ COMPLETED    module load fastp/0.20.1
#+ COMPLETED    
#+ COMPLETED    echo "Trimming downloaded datasets with fastp."
#+ COMPLETED    echo ""
#+ COMPLETED    
#+ COMPLETED    function trim_reads {
#+ COMPLETED    cd ${d}/$1/raw_reads
#+ COMPLETED    ls *.fastq.gz | cut -d "." -f "1" | sort > fastq_list
#+ COMPLETED    while read sample; do fastp -w ${t} -i ${sample}.fastq.gz -o ../cleaned_reads/${sample}_cleaned.fastq.gz
#+ COMPLETED    done<fastq_list
#+ COMPLETED    cd ..
#+ COMPLETED    }
#+ COMPLETED    
#+ COMPLETED    trim_reads Illumina
#+ COMPLETED    trim_reads Nanopore
#+ COMPLETED    
#+ COMPLETED      
#+ COMPLETED    module unload fastp/0.20.1
    
######################
# Indexing Reference #
######################

module load bwa/2020_03_19

echo "Indexing Reference"
echo ""
cd ${d}
bwa index ${g}

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
while read sample; do bwa mem ${d}/${g} ${sample}_cleaned.fastq.gz > ${sample}_mapped.sam
samtools sort ${sample}_mapped.sam > ${sample}_sorted.bam -@ ${t}
done<cleaned_list
cd ..
}

# align_reads Illumina
align_reads Nanopore

module unload bwa/2020_03_19
module unload samtools/1.16

#+  echo "Aligning Illumina datasets againts reference with minimap, using n threads."
#+  echo ""
#+  cd Nanopore
#+  b= ls -1 *.fastq.gz.fastp.gz
#+  for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr ${g} ${b} > ${b}.sam -t ${t}
#+  samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${t}
#+  rm ${b}.sam
#+  rm ${b}
#+  done
#+  cd ..
#+  
#+  module unload minimap2/2.24
#+  
#+  ###################################
#+  # Renaming and Indexing BAM files #
#+  ###################################
#+  
#+  module load samtools/1.16
#+  
#+  echo "Renaming files in bash"
#+  echo ""
#+  cd Illumina
#+  for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done
#+  cd ..
#+  f= ls -1 *.bam
#+  for f in *.bam; do samtools index ${f}; done
#+  cd Nanopore
#+  for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done
#+  cd ..
#+  f= ls -1 *.bam
#+  for f in *.bam; do samtools index ${f}; done
#+  
#+  module unload samtools/1.16
#+  
#+  
#+  ########## Done ###########
#+  
#+  end=`date +%s`
#+  elapsed=`expr $end - $begin`
#+  echo Time taken: $elapsed
#+  
#+  mdate=`date +'%d/%m/%Y %H:%M:%S'`
#+  # NOTE: If you are running a mac and having trouble with the code below,
#+  # ----- try using 'vm_stat' instead of 'vmstat'
#+  mcpu=$[100-$(vmstat 1 2|tail -1|awk '{print $15}')]%
#+  mmem=`free | grep Mem | awk '{print $3/$2 * 100.0}'`
#+  echo "$mdate | $mcpu | $mmem" >> ./stats-cpu
#+  ###############################################################
#+  #
} | tee logfile
