# mapping_cleaned_reads
Orientation to reference-based short-read mapping for advanced bioinformatics course at Utah Tech University

---

# Contents

-   [Objectives](#objectives)
-   [Genomic Filetypes](#genomic-filetypes)
-   [Basic Processing Steps](#basic-processing-steps)
-   [Exercise](#exercise)

---

# <a name="objectives"></a>
# Objectives

-  Understand why read mapping is an important part of bioinformatics pipelines
-  Understand the file conversions
---

# <a name="getting-set-up"></a>
# Getting set up
If you are here as a UTU student taking BIOL 4300, you should do the following:

1.  Login to your [Github](https://github.com/) account.

1.  Fork [this repository](https://github.com/KLab-UT/mapping-cleaned-reads), by
    clicking the 'Fork' button on the upper right of the page.

    After a few seconds, you should be looking at *your*
    copy of the repo in your own Github account.

1.  Click the 'Clone or download' button, and copy the URL of the repo via the
    'copy to clipboard' button. **note: if you have an SSH key with your github account, make sure you select the ```SSH``` tab**

1.  Login to the lonepeak cluster (CHPC) from your terminal

1.  On your lonepeak login node, navigate to where you want to keep this repo (
    I recommend having an exercise folder where you can clone repositories for the
    coding exercises). Then type:

        $ git clone the-url-you-just-copied

    and hit enter to clone the repository. Make sure you are cloning **your**
    fork of this repo.

1.  Next, `cd` into the directory:

        $ cd the-name-of-directory-you-just-cloned

1.  At this point, you should be in your own local copy of the repository.

    As you work on the exercise below, be sure to frequently `commit` your work
    and `push` changes to the *remote* copy of the repo hosted on Github.
---

# <a name="what-is-trimming"></a>
# What is mapping?

Doing genomic analyses is comparable to a 10 yr old child who accidentally shattered their mom's priceless artisan vase and is now frantically trying to put it back together. Only rather than having a few hundred pieces to deal with, a genomicists frequently work with millions of pieces. This is made easier when they have a reference to help orient the pieces (e.g., imagine the previously referenced 10-yr old found a photo of the vase in its original condition). In this tutorial, you will learn how to use bioinformatics tools to "map" reads (the broken pieces) to a reference genome (the vase photo).

Read files are usually stored within the ```fastq``` filetype. Before mapping, these reads should have been "cleaned" to remove adapter sequences and filter low-quality reads/bases (for an introduction to read cleaning, see the [trimming_raw_reads](https://github.com/KLab-UT/trimming_raw_reads) repository). A mapping (a.k.a. "aligning") software package can then be used to orient the reads to a reference. Ideally this reference is a well-annotated genome, but sometimes that is not available. In these scenarios, a bioinformatician may need to use an unannotated genome, the genome of a closely-related organism, or attempt to perform de-novo assembly to create a pseudo-reference using the read files. Mapped reads are stored in the ```.sam``` (or compressed vesion, ```.bam```) format (for an overview of genomic filetypes, see the [genomics-pipeline-intro](https://github.com/KLab-UT/genomics-pipeline-intro)).

Mapping reads is based on the evolutionary principle of homology. Homologous loci share evolutionary history. When mapping reads to a reference genome, we look for homologous regions (loci that share characteristics due to common ancestry). However, mutations cause homologous loci to have different nucleotide sequences. These mutations can be base substitutions (e.g., SNPs) or inversions/deletions. Mutations can lead to mismatches (the reference homolog doesn't exactly match the read due to a substitution) and/or gaps (a deletion event would result in a gap within the read that is not in the reference). Mapping software packages have to deal with these mismatches.

In addition to mismatches, data that comes from mRNA (e.g., RNA-seq) will contain reads that are not found in the genome due to splicing. In other words, genes contain coding regions (exons) which are interrupted by non-coding regions (introns). In the mRNA, these exons are joined together in a single genetic sequence. However, if this mRNA is being mapped to a reference genome, then accounting for the introns is necessary for accurate mapping. **When working with RNA-seq data, it is important to use a "splice-aware" aligner when mapping reads to a reference genome. Examples of splice-aware aligners include [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/), and [minimap2](https://github.com/lh3/minimap2).**

We will be using the software package ***[BWA](https://bio-bwa.sourceforge.net/bwa.shtml)***.

# Looking at the alignment program

The ```bwa``` help menu on the command line isn't great, so to see options for running the software see the [website](https://bio-bwa.sourceforge.net/bwa.shtml).

### Indexing
Before mapping reads to the reference, an "index" that breaks the reference into manageable chunks needs to be created. Similar to when you want to search for a term in a textbook, having a searchable index is much more efficient than looking through all of the material.

While in the ```Reference``` directory, create an index for the "GCF_000688415.1_ASM68841v1_genomic.fna" (the sars-cov-2 virus genome):

```
cd Reference
module load bwa
bwa index GCF_000688415.1_ASM68841v1_genomic.fna
```

You'll see that this creates a file called "covid-19-refseq.fasta.fai". This file contains 5 tab-delimited files. For organisms with multiple chromosomes, there will be a line for each chromosome

### Mapping

After creating the reference index, you are ready to map reads to the reference! We'll be using the ```mem``` algorithm. This version of ```bwa``` is efficient and accurate at mapping reads for long and short reads to large and small reference sequences for both paired-end and sing-end sequence data. The syntax is simple, try the following while in the main directory:

```
bwa mem Reference/GCF_000688415.1_ASM68841v1_genomic.fna example.fastq > example.sam
```

> To review the components of a ```.sam``` file, revisit the [genomics-processing-intro](https://github.com/KLab-UT/genomics-pipeline-intro) repository.

### Alignment File Processing

Alignment files in the ```.sam``` format can be very large. You can convert the ```example.sam``` file to a ```.bam``` (binary alignment map) using the software package [samtools](http://www.htslib.org/).

```
module load samtools
samtools sort example.sam > example.bam
```

Compare the file sizes between your ```.sam``` and ```.bam``` files:

```
ls example* -lh
```

You can also calculate your read depth at positions along your reference genome using samtools. Try the following:

```
samtools depth example.bam
```

The output of this command gives (1) the sequence name of the reference, (2) the position along the reference, and (3) the number of reads that mapped to that position. You'll notice from the output of this command that all positions have at least a depth of 1, and that many of the positions of the reference are left out. This is because ```depth``` only returns the positions that were mapped to. To see all of the positions along the reference genome, you can do the following command:

```
samtools depth -a example.bam
```

To get an average coverage per position in the reference sequence, you can perform the following command:

```
samtools depth  example.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

Once again, be careful! This only considers the sites along the reference with a minimum coverage of ```1``` when calculating the average. Can you think of a way to calculate average while incorporating all of the sites (including those with mapping depth of ```0```? (Hint: look at what we did before this section).

### Exercise
Run the script called ```mapping_cleaned_reads.sh``` (contained in the  bash_scripts directory). You can see flag and argument requirements by running the following command:

```
bash bash_scripts/mapping_cleaned_reads.sh -h
```

The path to the repository (```-d```) should be the path to your cloned
repository of your fork. The path to the working directory (```-w```) should be
the path to your personal directory in scratch. The illumina (```-I```) and
nanopore (```-N```) flags are for the text files that contain the references to
raw data. The reference genome (```-g```) is the fasta file contained within the
Reference repository. Your number of threads (```-t```) is the number of CPU
processors that you define in the batch script that runs your job.

> NOTE: You need to create a personal directory in scratch at this location:

```
/scratch/general/nfs1/
```

> For example, if my username was "u0123456", then I would create a directory in
this path called "u0123456". My path would look like this:

```
/scratch/general/nfs1/u0123456
```

Here is an example of how to run the script from the main directory of the repository using all of the flags:

```
bash bash_scripts/mapping_cleaned_reads.sh -I HumanNasalMicrobiota_Illumina.txt -N HumanNasalMicrobiota_Nanopore.txt -g GCF_000688415.1_ASM68841v1_genomic.fna -d ~/BIOL_4310/Exercises/Exercise_5/mapping_cleaned_reads -w /scratch/general/nfs1/u0123456 -t 4
```

Create a batch script to run the bash script with appropriate parameters for each flag. Then complete the worksheet before you push.

> Tip: You should do a trial run on an abbreviated dataset during an interactivesession. Create modified versions of  the Illumina and Nanopore text files so theycontain only one ID. Initiate an interactive session by running ```salloc```.Then you can execute the above bash command without submitting a batch job.

```
add worksheet.md logfile
git commit -m "ran script and answered worksheet questions"
git push
```

