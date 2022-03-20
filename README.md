# Combinatorial Indexing Analysis Pipeline

Pipeline for demultiplexing Combinatorial Indexing data of single cell Epitope measurements and generating a CELL * PROTEIN matrix.
Any arbitrary barcode pattern can be mapped to the reads, where the pattern can include: 
  - CI-barcodes
  - any Linker sequenes between barcodes
  - UMI sequence
  - gRNA sequence
  ( - RNA sequence coming soon ...)

The barcoding pattern is handed to the tool by a regex-like input parameter which summarizes the pattern sequence. E.g. [NNNNN][XXXXXXXXX][AGCTCATCGAC] is a barcoding pattern that contains three sequences: a CI-barcode [N...] (where the possibilities must be listed in an additional parameter), UMI sequence [X...] and a constant Linker sequence.

The Pipeline can be used with a different number of allowed mismatches in each sequence of the pattern. Each sequence of the pattern is mapped sequentially to the fastq-reads using semi-global alignments with Levenshtein distance.
Finally a UMI correction step can be performed to count UMIs within close proximity once.
(A version for counting RNA instead of Epitopes is in development).

A short overview of the Pipeline:
![Pipeline](https://github.com/tstohn/CombinatorialIndexingPipeline/blob/CITool/docs/media/PipelineReview.png)

The repository contains a few cpp tools that can be used for demultiplexing/ protein number counting seperately.
Otherwise you can also run the whole pipeline as a php script, which will perform demultiplexing & subsequent read counting.

Input are the raw fastq(.gz) files (the pipeline supports single or paired-end reads). However a single read is recommended if you want to run the Pipeline with a predefined number of mismatches in the overlapping region (stitch e.g. with fastq-join). In single read mode one pattern sequence after the other is sequentially mapped to the reads (e.g. first UMI pattern, then AB pattern as in image above), however in paired-end mode it might well be that the pattern sequence in the middle can not be completely mapped in ether read (forward & reverse), in that scenario this sequence is skipped as long as it is only a Linker sequence and ehter way not of interest for the CI Analysis. This however means that we can not assure the maximum number of mismatches in this region that the tool considers.
Output is a tsv file, with a column for the [protein], the [single cell ID], the [protein count] and dependant on the input parameters also a treatment of this cell and/or the cell origin (gRNA).


# Overview Pipeline/ Tools:

After compilation of the tools you can run *php src/Pipelines/analyze.php* to run the pipeline.

After compilation tools are found in *./bin*
  - **Demultiplexing**: Splitting the fastq-reads into tab seperated sequences. in the order of the barcode pattern
  - **READ PROCESSING**: Generating a Cell * Gene Matrix for the mapped reads
   
  
# Get started:

Firstly install necessary dependancies and build tools:
```bash
  make install
  make demultiplexing
  make processing
```
  
 For help how to run the pipeline run:
 ```bash
 php src/Pipelines/analyze.php --help
 ```

# Example

The Barcode File might look like this:

```
AGCTTAGC,ACGTTAAT,ATGCATGC
ACGTTAGC,ACTGCGAT,ACTGGATA
ACGT,ACTG
```
