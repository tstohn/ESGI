# ESGI
# Efficient Splitting of Generic Indices

[![LINUX BUILD (22.04)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/makefile.yml)
[![MACOS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/macOS.yml)
[![WINDOWS BUILD](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml/badge.svg)](https://github.com/tstohn/CombinatorialIndexingPipeline/actions/workflows/windows.yml)


<img src=https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/SYMBOL.png width="200" />

Pipeline for demultiplexing and counting generic barcoded sequencing data. Examples of technologies that can be demultiplexed are SPLiT-seq, Phospho-seq, SIGNAL-seq, scID-seq, spatial sequencing data and many,many more single-cell sequencing technologies. After demultiplexing ESGI runs UMI collapsing and creates a single-cell * feature matrix and can count any barcoded modality like RNA, antibody-tagged sequences, etc. 
Any arbitrary barcode pattern can be mapped to the reads, where the pattern can include: 
  - **Barcode elements**, e.g., for combinatorial-indexing like in split-and-pool approaches like SPLiT-seq, where a combination of different barcodes define a single cell. This barcode can also be a feature like in antibody-tagged sequences like in CITE-seq
  - **Constant linker elements**, e.g. linker sequences between barcodes of a constant known nucleotide-sequence
  - **UMIs**, or also several UMI barcodes within one read

ESGI can demultiplex also sequences where barcodes are of variable length (like staggers, where for a barcode at a specific position barcodes of different length are possible). ESGI can handle insertions,deletions and substitutions making it possible to demultiplex also erroneous data and *demultiplex* can handle several possible barcode-sequences in the same experiment: if one fastq-file contains different modalities with different barcode patterns. *demultiplex* then compares every fastq-line against all possible patterns and retains the best one, if it is uniquely the best and within the allowed mismatches.

The barcoding pattern is handed to the tool by a regex-like input parameter which summarizes the pattern sequence. E.g. [BC1.txt][10X][AGCTCATCGAC][BC2.txt] is a barcoding pattern that contains three sequences: a barcode element where all possible sequences are listed in the file BC1.txt (comma separated barcodes, they can be of different length), an UMI sequence of 10 random barcodes, a constant barcode with the sequence AGCTCATCGAC and finally a last barcode element with barcodes from the list in BC2.txt.

The Pipeline allows to set different mismatches for every barcode in the pattern: imagine BC2 has many more mismatches for some reason, or has a longer sequence and we can allow for more errors. We can then set the mismatches for BC2 higher than for other barcodes.

A short overview of ESGI:
![Pipeline](https://github.com/tstohn/CombinatorialIndexingPipeline/blob/master/docs/media/Overview.png)


# Overview ESGI:

__ESGI__   essentially performs two tasks: demultiplexing and counting. The tool *demultiplex* demultiplexes the input FASTQ-reads, essentially 'cutting' them into the dedicated barcode-sequences and the tool *count* creates a single-cell * feature matrix by counting the features like antibody-barcodes/RNA-sequences with UMI-collapsing. If barcode-sequences contain a RNA sequence ESGI can call STAR to align this part to a reference genome and then annotates the output of *demultiplex* with the STAR-result before running *count*. Next to *ESGI* this repo provides all three tools: *demultiplex*, *count* and *annotate* (although annotate is only used to annotate the output of demultiplex with the mapped STAR reads and probably of little interest to users). *ESGI* can process fastq(.gz) files, supports single or paired-end reads or also simple txt files with nucleotide-sequences in every line.

__Input__   is an ini-file (like [example.ini](https://github.com/tstohn/ESGI/blob/master/src/test/test_data/test_esgi/esgi_example.ini)) which contains all the necessary information about the input-fastqs, the barcode-patterns, which barcodes define a single-cell, allowed mismatches, etc. See below section 'Example usage', run the tools with the --help flag or visit our website for more information.

__Output__   is a tsv file, with a column for the *feature ID*, the *single-cell ID*, and the UMI-collpased *count*. If additional annotation information was given (like barcodes that encode batches, treatment-conditions, spatial-coordinates, etc. this information is present in additional columns). Additionally, ESGI provides UMI-information (amplification of different UMIs), information on errors (insertions, deletions, substitutions) for the different barcode positions, the simple demultiplexed fastq-line output (the fastq-sequences 'cut' into its barcode-pattern) and other useful information.

ESGI can be run as a wrapper of *demultiplex* and *count* with an ini-file that contains all the information or its tools can be run individually by providing the input-information in the command line.

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

# Example usage:

You can run ESGI with an initialization-file (file-ending=.ini) that contains all the information about the experiment like:
```bash
  ./bin/esgi myExperiment.ini
```

The myExperiment.ini could look like this:

Next to the input files the ini needs to state at least:
- the pattern-file: a txt file that contains the pattern(s) that should be matched to the input sequences 
- the mismatch-file: a txt file with comma-seperated mismatched for every element in the pattern-file
- the IDs for single-cells, features, UMIs if present (these are the positions of elements in the pattern-file, where indexing starts from 0)
  
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

Below are examples of the additional files given in the ini:

The pattern that should be mapped to the sequences is written into the pattern.txt file.
The file can contain one or several patterns, every additional pattern must be written in a new row.
This file lists all pattern-elements like barcodes, UMIs, etc. and those elements are enclosed by squared brackets.
The pattern can have a name (e.g., PATTERN_NAME:). This is not required, but can be handy if several patterns are provided, as ESGI creates one file of demultiplexed reads for every pattern.
Possible elements inlude:
  - a constant nucleotide sequence given as string of A,G,C,Ts, e.g. [GCATTACG].
  - barcode sequences, that are given by the path to a txt file. This file contains a comma-seperated list of all possible barcodes at this position.
  - UMI stated as <number>X, e.g. [15X]. You can also use this if you do not care about constant sequences. E.g., imagine we have the pattern [barcodes.txt][AAAA][barcodes.txt] but we do not care about the [AAAA] at all, we can simply use [barcodes.txt][4X][barcodes.txt]. This pattern element [4X] does not have to be used as UMI. When running ESGI or count we state the index in the pattern that we want to be used as UMI, this pattern element has to be a <number>X, but not every <number>X must be used as UMI. You can also ahve several <number>X elements and use them all as UMI, then ESGI/ count concatenates all the <number>X elements that are used as UMI and uses them as one long UMI.
  - genomic sequences like RNA/DNA, that need to be aligned to a reference genome with STAR, are listed as [DNA].
ESGI makes use of two additional elements for special barcoding cases:
  - [-] seperates forward/reverse read strictly. The pattern generally covers the forward and reverse read (assuming reverse complements of the reverse read).
    Sometimes the user might not have an overlap and wants to strictly seperate the pattern for the forward and reverse read. In that case [-] can be used. Be aware
    that ESGI expects per default an overlap of forward and reverse read, therefore, even when supplying the [-] element, the elements after [-] would map to the end of the           reverse complement of the reverse read. If this is not wanted the user can set the indepent-parameter (run esgi --help for more info).
    Example: [ACCAGT][-][BC_reverse.txt][AAA]. We would expect in the reverse-reads following sequence: TTT followed by the reverse complement of a sequence in BC_reverse.
    If we add a line independent=1 in the ini we expect the reverse-read to start with a barcode from BC_reverse (NOT the reverse complement of the barcode) followed by AAA.
 - [*] says we simply stop mapping when encountering this element from the forward or reverse read. Can be handy if there is a sequqence in the middle of the read that we are not    interested in.

```txt
PATTERN_NAME:[GCATTACG][/USER/DATA/MYEXPERIMENT/BC1.txt][CAGTACCG][/USER/DATA/MYEXPERIMENT/ANTIBODY_BC.txt][10X][/USER/DATA/MYEXPERIMENT/BC2.txt]
```
mismatches.txt
This file states the allowed number of mismatches in every pattern and for the UMI with how many mismatches UMIs are corrected.
(2 MM in constant elements, 1MM in barcode elements and aligning UMIs with 1MM)
```txt
2,1,2,1,1,1
```

Finally, some examples for the files containing the barcodes, or names for features, conditions, that can be additionally assigned if needed.
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
# Running ESGI with staggers

ESGI can demultiplex patterns with staggers (barcodes of variable length). One example pattern would be [A|AC|ACG|ACGT][GGGG] where we expect first a barcode of length 1 to 4 followed by a constant element GGGG. The problem of staggers is that several barcodes might map equally well: imagine we have the read ACGGGG. Now we first map the stagger barcode and A,AC,ACG all would map equally well (when mapping barcodes of different length ESGI does not punish deletions at the end of the barcode). Therefore, we would ether have to map the whole pattern first to see that actually AC and GGGG would be the best split, or at least map the stagger sequence together with the barcode that follows the stagger!!

ESGI has two features that makes it possible to match staggers: (1) Barcodes in barcode-elements (elements described by a txt file that contains all possible barcodes) can have variable length and (2) ESGI can demultiplex several patterns simultaneously.

1.) use a single pattern and merge the stagger with the constant sequence. If we would not merge them and have a pattern with the stagger and the constant element many reads would be discarded because of ambiguous mapping, since if a read contains barcode 'ACG' also barcode 'A' and 'AC' would map and ESGI would discard the read (ESGI does not look ahead, but at every barcode position tries to find the best match, if there are several matches the read is discarded). But we can merge the stagger with the constant sequence and even allow for a mismatch, ESGI would still find the pattern that at least matches best.

pattern.txt
```txt
[stagger_barcodes.txt]
```
stagger_barcodes.txt
```txt
AGGGG,ACGGGG,ACGGGGG,ACGTGGGG
```

2.) The second option would be to describe an individual pattern for every stagger, and allow only for no or very little mismatches in the staggers with very few nucleotides. This way we prevent to map a wrong barcode with insertions/deletions to a stagger. Additionally, you could map with hamming distance only in the barcodes with the -H flag.

pattern.txt
```txt
PATTERN_1:[A][GGGG]
PATTERN_2:[AC][GGGG]
PATTERN_3:[ACG][GGGG]
PATTERN_4:[ACGT][GGGG]
```
mismatches.txt
```txt
0,1
0,1
1,1
1,1
```


# Points to consider

- at the moment we do not compile ESGI with htslib under Windows, you need to build it yourself if you want to run ESGI with STAR (and annotate) on Windows. You can run ESGI without problem on Windows without RNA-mapping - or even run only demultiplex with an RNA-sequence if you just want to split the barcode-sequences. If you install htslib on Windows and want to compile it with htslib, set the variable HTSlib_AVAILABLE:=yes in the Makefile.
- If you download the binaries you can run ESGI right-away without RNA-mapping. If you want to run RNA-mapping (with STAR and annotate) you need to manually install STAR and htslib. To check if it works you can run
```bash
  make test_esgi_RNA
```
- at the moment the multi-pattern option is supported by ESGI only if all patterns belong to the same modality (this is because *count* is called only once in ESGI for all patterns together) and the number of columns, positions of single-cell barcodes, etc. is the same in all patterns. An example for this could be if the data contains reads with different barcodes at the barcode elements positions, different barcodes have different constant linkers attached, or certain barcodes at one positions were only combined with certain other barcodes at a later position etc. and the user wants to explicitely state these different patterns (for more details see website). Still all barcodes are at the same position and can all be counted together. If you have a FASTQ with several modalities/ a hierachical pattern we recommend to run ESGIs tools individually: 1.) run *demultiplex* in multi-pattern mode (see demultiplex --help) and for the individual outputs of *demultiplex* (one for every pattern) run *count*.


