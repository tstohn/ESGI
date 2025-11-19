# ESGI
# Efficient Splitting of Generic Indices

[![LINUX BUILD (22.04)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml)
[![MACOS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml)
[![WINDOWS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml)


<img src=https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/DemultiplexingLogo.png width="200" />

Pipeline for demultiplexing and counting generic barcoded sequencing data. Examples of technologies that can be demultiplexed are SPLiT-seq, Phospho-seq, SIGNAL-seq, scID-seq, spatial sequencing data and many,many more single-cell sequencing technologies. After demultiplexing ESGI runs UMI collapsing and creates a single-cell * feature matrix and can count any barcoded modality like RNA, antibody-tagged sequences, etc. 
Any arbitrary barcode pattern can be mapped to the reads, where the pattern can include: 
  - **Variable Barcodes**, e.g., for combinatorial-indexing like in split-and-pool approaches like SPLiT-seq, where a combination of variable barcodes define a single cell. This barcode can also be a feature like in antibody-tagged sequences like in CITE-seq
  - **Constant Barcodes**, e.g. linker sequences between barcodes of a non-variable known nucleotide-sequence
  - **UMI Barcodes**, or also several UMI barcodes within one read

ESGI can demultiplex also sequences where barcodes are of variable length (like staggers, where for a barcode at a specific position barcodes of different length are possible). ESGI can handle insertions,deletions and substitutions making it possible to demultiplex also erroneous data and *demultiplex* handle several possible barcode-sequences in the same experiment: if one fastq-file contains different modalities with different barcode patterns the ESGI can take several barcode-sequence patterns as input searching for them simultaneously in the fastq.

The barcoding pattern is handed to the tool by a regex-like input parameter which summarizes the pattern sequence. E.g. [BC1.txt][10X][AGCTCATCGAC][BC2.txt] is a barcoding pattern that contains three sequences: a variable-barcode where all possible sequences are listed in the file BC1.txt (comma separated barcodes, they can be of different length), an UMI sequence of 10 random barcodes, a constant barcode with the sequence AGCTCATCGAC and finally a last variable barcode from the list in BC2.txt.

The Pipeline allows to set different mismatches for every barcode in the pattern: imagine BC2 has many more mismatches for some reason, or has a longer sequence and we can allow for more errors. We can then set the mismatches for BC2 higher than for other barcodes.

A short overview of ESGI:
![Pipeline](https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/Figure_1.png)


# Overview ESGI:

__ESGI__   essentially performs two tasks: demultiplexing and counting. The tool *demultiplex* demultiplexes the input fastq-reads, essentially 'cutting' them into the dedicated barcode-sequences and the tool *count* creates a single-cell * feature matrix by counting the features like antibody-barcodes/RNA-sequences with UMI-collapsing. If barcode-sequences contain a RNA sequence ESGI can call STAR to align this part to a reference genome and then annotates the output of *demultiplex* with the STAR-result before running *count*. Next to *ESGI* this repo provides all three tools: *demultiplex*, *count* and *annotate* (although annotate is only used to annotate the output of demultiplex with the mapped STAR reads and probably of little interest to users). *ESGI* can process fastq(.gz) files, supports single or paired-end reads or also simple txt files with nucleotide-sequences in every line.

__Input__   is an ini-file (like [example.ini](https://github.com/tstohn/ESGI/blob/master/src/test/test_data/test_esgi/esgi_example.ini)) which contains all the necessary information about the input-fastqs, the barcode-patterns, which barcodes define a single-cell, allowed mismatches, etc. See below section 'Example usage', run the tools with the --help flag or visit our website for more information.

__Output__   is a tsv file, with a column for the *feature ID*, the *single-cell ID*, and the UMI-collpased *count*. If additional annotation information was given (like barcodes that encode batches, treatment-conditions, spatial-coordinates, etc. this information is present in additional columns). Additionally, ESGI provides UMI-information (amplification of different UMIs), information on errors (insertions, deletions, substitutions) for the different barcode positions, the simple demultiplexed fastq-line output (the fastq-sequences 'cut' into its barcode-pattern) and other useful information.

# Documentation

Check out the [ESGI website](https://tstohn.github.io/ESGI.documentation/) for more detailed information with examples and use-cases how to run ESGI and its tools demultiplex/count.
You can also run:
```bash
  esgi --help
```
or for any of the other tools:
```bash
  demultiplex --help
  count --help
```

It might also help to look into the Makefile, browse through some tests there or look at the analyses that were run for the paper in [ESGI paper analyses]([https://github.com/username/repository)(https://github.com/tstohn/Analysis-EZGI))
  
# Getting started:

The easiest way to get started right away is to simply download the binaries of ESGI for your system (Windows, Linux or Mac) here [RELEASE](https://github.com/tstohn/ESGI/releases).
They contain the whole repository. __ESGI__, __demultiplex__, __annotate__ and __count__ are in bin and you can find the library for ESGI in lib if you want to develop with it yourself. 
The releases do not contain htslib - which is needed for *annotate* if you want to map RNA-sequences - and STAR. Therefore, both have to be downloaded manually. If your sequences contain no RNA sequences you can work with the release right away.

Otherwise, you can also build ESGI on your own if you prefer (or if integrating htslib into the release for RNA-mapping does not work).
Therefore, you can install necessary dependencies and then build ESGI (and its tools):
```bash
  make install
  make esgi
```

The individual tools can also be build separately:
```bash
  make demultiplex
  make count
```

If you prefer to not run ESGI inside the terminal you can download the ESGI-APP with a graphical user interface here: [APP](https://github.com/tstohn/ESGI-GUI). It contains an interface to run demultiplex and count on your data. The app is only tested on MacOS, is intended to run on any standard laptop and does not provide STAR-mapping of RNA-sequencing.

# Example usage:

You can run ESGI with an initialization-file (file-ending=.ini) that contains all the information about the experiment like:
```bash
  ./bin/esgi myExperiment.ini
```

The myExperiment.ini could look like this:
```ini
  # THIS IS A MINIMAL EXAMPLE FOR THE .INI FILE
  # THE FILE IS NOT WHITESPACE SENSITIVE
  # FOR MORE DETAILS RUN: ESGI --help

  # INPUT FILES: can be fastq or txt and raw or gzipped
  forward=/USER/DATA/MYEXPERIMENT/INPUT/forward_reads.fastq.gz
  #reverse file is optionally in case we have fw and rv reads
  reverse=/USER/DATA/MYEXPERIMENT/INPUT/reverse_reads.fastq.gz

  #output directory
  output=/USER/DATA/MYEXPERIMENT/OUTPUT

  pattern=/USER/DATA/MYEXPERIMENT/pattern.txt
  mismatches=/USER/DATA/MYEXPERIMENT/mismatches.txt

  # Barcodes that are used to define individual-cells (can be one barcode or a combination for,e.g. split-and-pool experiments)
  # the indices for single cells are the barcode-positions in the pattern.txt file
  # indexing starts at 0 and counts EVERY pattern element that is defined by '[]', even [-] or [*]
  SC_ID=1,5
  FEATURE_ID=3
  # feature name is optionally if we have, e.g. antibody-barcodes that should be assinged names like AGCAGCAT-> antibody of EGFR
  FEATURE_NAMES=/USER/DATA/MYEXPERIMENT/FEATURES.txt
  # barcode-file that is also used for additional annotations, e.g. if cells that were in certain wells during indexing in barcode-round 1 (BC1.txt) were treated differently
  # (this can also be a list (comma or whitespace separated) if we have several annotations - like another barcode encoding something else)
  ANNOTATION_IDs=/USER/DATA/MYEXPERIMENT/BC1.txt
  # file containing names for the annotation-barcodes above like treatment-condition names
  # (this can also be a list (comma or whitespace separated) if we have several annotations - like another barcode encoding something else)
  ANNOTATION_NAMES=/USER/DATA/MYEXPERIMENT/TREATMENTS.txt
  UMI_ID=4
  
  threads=10
  prefix=MYEXPERIMENT

```

with additional files looking like this:

pattern.txt
(variable-barcodes are defined in .txt files, UMIs are defined by [<base-number>X] - for more detail look into the documentation)
```txt
PATTERN_NAME:[GCATTACG][/USER/DATA/MYEXPERIMENT/BC1.txt][CAGTACCG][/USER/DATA/MYEXPERIMENT/ANTIBODY_BC.txt][10X][/USER/DATA/MYEXPERIMENT/BC2.txt]
```
mismatches.txt
(2 MM in constant-barcodes, 1MM in variable-barcodes and aligning UMIs with 1MM)
```txt
2,1,2,1,1,1
```

BC1.txt
```txt
AC,CACA,GACTGA,GAACTGAA
```
BC2.txt
```txt
ATAT,CGAT,TAAG,CCGG
```
ANTIBODY_BC.txt
```txt
AAAA,CCCC,TTTT,GGGG
```

FEATURES.txt
```txt
pEGFR,pRAS,pMEK,pERK
```
TREATMENTS.txt
```txt
CONTROL,EGFRi,CONTROL,EGFRi
```

# Points to consider

- at the moment we do not compile ESGI with htslib under Windows, you need to build it yourself if you want to run ESGI with STAR (and annotate) on Windows. You can run ESGI without problem on Windows without RNA-mapping - or even run only demultiplex with an RNA-sequence if you just want to split the barcode-sequences. If you install htslib on Windows and want to compile it with htslib, set the variable HTSlib_AVAILABLE:=yes in the Makefile.
- If you download the binaries you can run ESGI right-away without RNA-mapping. If you want to run RNA-mapping (with STAR and annotate) you need to manually install STAR and htslib. To check if it works you can run
```bash
  make test_esgi_RNA
```
- at the moment ESGI does not support the multi-pattern option. If you have a fastq with many different modalities we recommend to run ESGIs tools individually: 1.) run *demultiplex* in multi-pattern mode (see demultiplex --help) and for the different outputs of demultiplex (one for every pattern) run *count*


