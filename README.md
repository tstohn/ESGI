# Combinatorial Indexing Pipeline

Pipeline for demultiplexing Combinatorial Indexing data of single cell Epitope measurements and generating a CELL * PROTEIN matrix.
Any arbitrary barcode pattern can be mapped to the reads, where the pattern can include: 
  - CI-barcodes
  - any Linker sequenes between barcodes
  - UMI sequence
  - gRNA sequence

The Pipeline can be used with a different number of allowed mismatches in each sequence of the pattern. Each sequence of the pattern is mapped sequentially to the fastq-reads using semi-global alignments with Levenshtein distance.
Finally a UMI correction step can be performed to count UMIs within close proximity once.

A short overview of the Pipeline:
![Pipeline](https://github.com/tstohn/CombinatorialIndexingPipeline/blob/CITool/docs/media/PipelineReview.png)

The repository contains a few cpp tools that can be used for demultiplexing/ protein number counting seperately.
Otherwise you can also run the whole pipeline as a php script, which will perform demultiplexing & subsequent read counting.


# Overview Pipeline/ Tools:

After compilation of the tools you can run *php src/Pipelines/analyze.php* to run the pipeline.

After compilation tools are found in *./bin*
  - **Demultiplexing**: Splitting the fastq-reads into tab seperated sequences. in the order of the barcode pattern
  - **READ PROCESSING**: Generating a Cell * Gene Matrix for the mapped reads
   
  
# Get started:

  `make install`  
  `make demultiplexing`  
  `make processing` 
  
 - run *php src/Pipelines/analyze.php --help* to get started
