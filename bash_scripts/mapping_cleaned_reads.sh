#!/bin/bash

{
usage="$(basename "$0") [-h] [-lI <SRA_list_Illumina>] [-lN <SRA_list_Nanopore>] [-t <threads>]
This program will call variants using freebayes in given SRA NGS sequences files to obtain major viral variants.
    -h  show this help text
    -lI  File or path to SRA accession list for Illumina data in tabular format
    -lN  File or path to SRA accession list for Nanopore data in tabular format
    -g  PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name.
    -t  Number of CPU processors"
options=':hl:g:a:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    lI) lI=$OPTARG;;
    lN) lN=$OPTARG;;
    g) g=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$lI" ] || [ ! "$lN" ] || [ ! "$g" ] || [ ! "$t" ]; then
  echo "arguments -lI, -lN, -g, and -t must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`

##################################################################################
#  Download fastq files for 22 samples from the Heikema et. al (2020) paper      # 
##################################################################################

module load sra-toolkit/3.0.2

mkdir -p Illumina
echo "Downloading Illumina SRA files from the given list of accessions"
cd Illumina
prefetch --max-size 800G -O ./ --option-file ${lI}
echo "Converting Illumina SRA files to fastq.gz"
ls -p | grep ERR > sra_dirs
while read i; do mv "$i"*.sra .; done<sra_dirs
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done
cd ..

mkdir -p Nanopore
echo "Downloading Nanopore SRA files from the given list of accessions"
cd Nanopore
prefetch --max-size 800G -O ./ --option-file ${lN}
echo "Converting Nanopore SRA files to fastq.gz"
ls -p | grep ERR > sra_dirs
while read i; do mv "$i"*.sra .; done<sra_dirs
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done
cd ..

module unload sra-toolkit/3.0.2

echo "SRA files were downloaded and .fastq files extracted in Illumina and Nanopore directories"
echo ""
echo "Done"
echo ""

##################################################################################
# Trimming downloaded Illumina datasets with fastp, using 16 threads (-w option) #
##################################################################################

module load fastp/0.20.1

echo "Trimming downloaded Illumina datasets with fastp."
echo ""

cd Illumina
z= ls -1 *.fastq.gz
for z in *.fastq.gz; do fastp -w ${t} -i ${z} -o ${z}.fastp
gzip ${z}.fastp
done
cd ..

cd Nanopore
z= ls -1 *.fastq.gz
for z in *.fastq.gz; do fastp -w ${t} -i ${z} -o ${z}.fastp
gzip ${z}.fastp
done
cd ..

module unload fastp/0.20.1

######################
# Indexing Reference #
######################

module load samtools/1.16

echo "Indexing Reference"
echo ""
samtools faidx ${g}

module unload samtools/1.16

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

module load minimap2/2.24

echo "Aligning Illumina datasets againts reference with minimap, using n threads."
echo ""
cd Illumina
b= ls -1 *.fastq.gz.fastp.gz
for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr ${g} ${b} > ${b}.sam -t ${t}
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${t}
rm ${b}.sam
rm ${b}
done
cd ..


echo "Aligning Illumina datasets againts reference with minimap, using n threads."
echo ""
cd Nanopore
b= ls -1 *.fastq.gz.fastp.gz
for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr ${g} ${b} > ${b}.sam -t ${t}
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${t}
rm ${b}.sam
rm ${b}
done
cd ..

module unload minimap2/2.24

###################################
# Renaming and Indexing BAM files #
###################################

module load samtools/1.16

echo "Renaming files in bash"
echo ""
cd Illumina
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done
cd ..
f= ls -1 *.bam
for f in *.bam; do samtools index ${f}; done
cd Nanopore
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done
cd ..
f= ls -1 *.bam
for f in *.bam; do samtools index ${f}; done

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
#
