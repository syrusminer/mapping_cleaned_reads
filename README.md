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

Read files are usually stored within the fastq filetype. Before mapping, these reads should have been "cleaned" to remove adapter sequences and filter low-quality reads/bases (for an introduction to read cleaning, see the [trimming_raw_reads](https://github.com/KLab-UT/trimming_raw_reads) repository). A mapping (a.k.a. "aligning") software package can then be used to orient the reads to a reference. Ideally this reference is a well-annotated genome, but sometimes that is not available. In these scenarios, a bioinformatician may need to use an unannotated genome, the genome of a closely-related organism, or attempt to perform de-novo assembly to create a pseudo-reference using the read files.

Mapping reads is based on the evolutionary principle of homology. Homologous loci share evolutionary history. When mapping reads to a reference genome, we look for homologous regions (loci that share characteristics due to common ancestry). However, mutations cause homologous loci to have different nucleotide sequences. These mutations can be base substitutions (e.g., SNPs) or inversions/deletions. Mutations can lead to mismatches (the reference homolog doesn't exactly match the read due to a substitution) and/or gaps (a deletion event would result in a gap within the read that is not in the reference). Mapping software packages have to deal with these mismatches.

In addition to mismatches, data that comes from mRNA (e.g., RNA-seq) will contain reads that are not found in the genome due to splicing. In other words, genes contain coding regions (exons) which are interrupted by non-coding regions (introns). In the mRNA, these exons are joined together in a single genetic sequence. However, if this mRNA is being mapped to a reference genome, then accounting for the introns is necessary for accurate mapping. **When working with RNA-seq data, it is important to use a "splice-aware" aligner when mapping reads to a reference genome. Examples of splice-aware aligners include [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/), and [minimap2](https://github.com/lh3/minimap2).**

Mapped reads are stored in the ```.sam``` (or compressed vesion, ```.bam```) format (for an overview of genomic filetypes, see the [genomics-pipeline-intro](https://github.com/KLab-UT/genomics-pipeline-intro)). We will be using the software package **"minimap2"**.

# Looking at an alignment program

While on the CHPC in this repository:
```
module load minimap2/2.24
minimap2 -h
```

From the minimap2 help menu, you will see sections on "Indexing", "Mapping", "Alignment", and "Input/Output". Let's break these down:

### Indexing
Before mapping reads to the reference, an "index" that breaks the reference into manageable chunks needs to be created. Similar to when you want to search for a term in a textbook, having a searchable index is much more efficient than looking through all of the material. When creating your index, you can specify how long each chunk will be (```-k```), 

Create an index for the 
```

add worksheet.md logfile
# Make sure you don't add all of the genomic files, just the worksheet.md
# The .gitignore includes the genomic files, so they won't be added by default
git commit -m "ran script and answered worksheet questions"
git push
```

