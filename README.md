# ESGI
# Efficient Splitting of Generic Indices

[![LINUX BUILD (22.04)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml)
[![MACOS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml)
[![WINDOWS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml)


<img src=https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/DemultiplexingLogo.png width="200" />

Pipeline for demultiplexing and counting generic barcoded sequencing data. Examples of technologies that can be demultiplexed are SPLiT-seq, Phospho-seq, SIGNAL-seq, scID-seq, spatial sequencing data and many,many more single-cell sequencing technologies. After demultiplexing ESGI runs UMI collapsing and creates a single-cell * feature matrix and can count any barcoded modality like RNA, antibody-tagged sequences, etc. 
Any arbitrary barcode pattern can be mapped to the reads, where the pattern can include: 
  - **Variable Barcodes**, e.g., for combinatorial-indexing like in split-and-pool approaches like SPLiT-seq where a combination of variable barcodes define a single cell. This barcode can also be a feature like in antibody-tagged sequences like in CITE-seq
  - **Constant Barcodes**, e.g. linker sequences between barcodes of a non-variable known nucleotide-sequence
  - **UMI Barcodes**, or also several UMI barcodes within one read

ESGI can demultiplex also sequences where barcodes are of variable length (like staggers, where for a barcode at a specific position barcodes of different length are possible). ESGI can handle insertions,deletions and substitutions making it possible to demultiplex also erroneous data and *demultiplex* handle several possible barcode-sequences in the same experiment: if one fastq-file contains different modalities with different barcode patterns the ESGI can take several barcode-sequence patterns as input searching for them simultaneously in the fastq.

The barcoding pattern is handed to the tool by a regex-like input parameter which summarizes the pattern sequence. E.g. [BC1.txt][10X][AGCTCATCGAC][BC2.txt] is a barcoding pattern that contains three sequences: a variable-barcode where all possible sequences are listed in the file BC1.txt (comma separated barcodes, they can be of different length), an UMI sequence of 10 random barcodes, a constant barcode with the sequence AGCTCATCGAC and finally a last variable barcode from the list in BC2.txt.

The Pipeline allows to set different mismatches for every barcode in the pattern: imagine BC2 has many more mismatches for some reason, or has a longer sequence and we can allow for more errors. We can then set the mismatches for BC2 higher than for other barcodes.

A short overview of ESGI:
![Pipeline](https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/ScDemultiplex.png)


# Overview ESGI:

__ESGI__ essentially performs two tasks: demultiplexing and counting. If barcode-sequences contain a DNA/RNA sequence ESGI can call STAR to align this part to a reference genome and then annotates the output of the demultiplexing with the STAR-result before running counting. Next to ESGI this repo provides all three tools: *demultiplex*, *count* and *annotate*. ESGI can process the raw fastq(.gz) files, supports single or paired-end reads or also simple txt files with nucleotide-sequences in every line.

__Input__ is an ini-file ([example.ini](https://github.com/tstohn/ESGI/blob/master/src/test/test_data/test_esgi/esgi_example.ini) which contains all the necessary information about the input-fastqs, the barcode-patterns, which barcodes define a single-cell, allowed mismatches, etc.
__Output__ is a tsv file, with a column for the *feature ID*, the *single-cell ID*, and the UMI-collpased *count*. Additionally ESGI provides UMI-information (amplification of different UMIs), information on errors (insertions, deletions, substitutions) for the different barcode positions, the simply demultiplexed fastq-line output (the fastq-sequences *cut* into its barcode-pattern) and other useful information.

# Documentation

Check out the [ESGI website](https://tstohn.github.io/ESGI.documentation/) for more detailed information with examples and use-cases how to run ESGI and its tools demultiplex/count.
You can also run:
```bash
  ESGI --help
```
or for any of the other tools:
```bash
  demultiplex --help
  count --help
  annotate --help
```

It might also help to look into the Makefile, browse through some tests there or look at the analyses that were run for the paper in [ESGI paper analyses]([https://github.com/username/repository)(https://github.com/tstohn/Analysis-EZGI))
  
# Getting started:

The easiest way to get started right away is to simply download the binaries of ESGI for your system (Windows, Linux or Mac) here [RELEASE](https://github.com/tstohn/ESGI/releases).
They contain the whole repository. __ESGI__, __demultiplex__, __annotate__ and __count__ are in bin and you can find the library for ESGI in lib if you want to develop with it yourself. 
The releases do not contain htslib - which is needed for *annotate* if you want to map dna/RNA-sequences - and STAR. Therefore, both have to be downloaded manually. If your sequences contain no DNA/RNA sequences you can work with the release right away.

Otherwise, you can also build ESGI on your own if you prefer (or if integrating htslib into the release for RNA-mapping does not work).
Therefore, you can install necessary dependencies and then build ESGI (and its tools):
```bash
  make install
  make ESGI
```

The individual tools can also be build and tested separately:
```bash
  make demultiplex
  make count
```

# Example usage:

You can run ESGI with an initialization-file (file-ending=.ini) that contains all the information about the experiment like:
```bash
  ./bin/ESGI myExperiment.ini
```

The myExperiment.ini could look like this:
```ini
  # THIS IS AN EXAMPLE FILE FOR THE INPUT.INI FILE FOR ESGI.
# THE FILE IS NOT WHITESPACE SENSITIVE
# SOME PARAMETERS ARE REQUIRED, OTHERS ARE NOT AND MARKED AS OPTIONAL

##############################
# MINIMAL ARGUMENTS
##############################

# INPUT FILES: can be fastq or txt and raw or gzipped
forward=src/test/test_data/test_esgi_RNA/SIGNALseq_with_toyRNA_1.fastq
#reverse file is optionally in case we have fw and rv reads
reverse=src/test/test_data/test_esgi_RNA/SIGNALseq_with_toyRNA_2.fastq

#output directory
output=bin

# a file containing all possible patterns that can exist, please make sure all files are in the same directory 
#(for compatability with the count tool which has only file names in the header of the inoput file and expects all barcodes in the same dir)
pattern=src/test/test_data/test_esgi_RNA/pattern_RNA.txt
# file containing all mismatches per pattern
mismatches=src/test/test_data/test_esgi_RNA/mismatches_RNA.txt

# PARAMETERS DEFINING SC_ID=Columns for single cells, FEATURE_ID for feature column, etc.
# indexing starts at 0 and taking EVERY pattern element in '[]' into consideration, even [-] or [*]
# at least SC_ID (barcode to distinguish individual cells) and FEATURE_ID must be given to count barcodes that
# encode a feature
SC_ID=2,4,6
FEATURE_ID=0
# feature name is optionally if we have, e.g. RNA or want simple barcodes as features without mapping a barcode to a name like AGCAGCAT->EGFR_ANITBODY
FEATURE_NAMES=
#ANNOTATION parameters are optionally if we have no additional annotations like treament conditions, batches, etc. encoded in a barcode
ANNOTATION_IDs=
#list of files, every file contains a list of annotations names in the same order as the barcodes for the specific annotation ID
ANNOTATION_NAMES=
#UMI_ID is optionally if we DO NOT want to collapse UMIs
UMI_ID=7

##############################
# ADDITIONAL ARGUMENTS
# boolean variables like 'independent' can be set with TRUE/FALSE or 0/1 or YES/NO
##############################

threads=1
# default 0, means the reverse read is a 'true' reverse and we sequentially map the pattern from 5' to 3' of forward read and then continue
# with reverse complement of reverse, set independent=1 if you want to continue mapping the reverse read independenlty from its own start (no reverse complement)
# in this case there has to be a [-] element in the pattern, everything before is mapped sequentially to forward, after sequentially to reverse
independent=0
#custom file prefix, e.g., date of run, name, etc.
#this prefix is added to the demultiplexed result, the tool count then add another prefix=COUNTDATA in front
prefix=SIGNALseq
writeFailedLines=1
writeStats=1
hamming=0
#default is input.threads * 100000
fastqReadBucketSize=500000
#filter to remove 'false' reads. It is a percentage x given as a value from [0-1[ 
#reads are only retained if >x percent of reads with the SAME UMI have the SAME FEATURE-BARCODE and SINGLE-CELL ID
umiThreshold=0.0
umiAbundance=0.0
umiCollapsing=1
SC_ID_string=0
# a file of barcodes that should be mapped to other barcodes, for more details look into the tool count
barcodeSharing=

##############################
# STAR ARGUMENTS
##############################
#optional path to the STAR executable to run
STAR=
#folder that holds the STAR genome index files, must be prebuild
genomeDir=src/test/test_data/test_esgi_RNA/star_index/toy
#the annotation of STAR that should be counted: GX=gene id, GN=gene name
feature=GX
```

# Points to consider

- at the moment we do not compile ESGI with htslib under Windows, you need to build it yourself if you want to run ESGI with STAR (and annotate) on Windows. You can run ESGI without problem on Windows without RNA-mapping - or even run only demultiplex with an RNA-seqeunce if you just want to split the barcode-sequences. If you install htslib on Windows and want to compile it with htslib, set the variable HTSlib_AVAILABLE:=yes in the Makefile.
- If you download the binaries you can run ESGI right-away without RNA-mapping. If you want to run RNA-mapping (with STAR and annotate) you need to manually install STAR and htslib. To check if it works you can run
```bash
  make test_esgi_RNA
```
- at the moment ESGI does not support the multi-pattern option. If you have a fastq with many different modalities we recommend to run ESGIs tools individually: 1.) run *demultiplex* in multi-pattern mode (see demultiplex --help) and for the different outputs of demultiplex (one for every pattern) run *count*


