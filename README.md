# SingleCellGenomeTools
Some simple tools to process fastq-files of nucleotide tagged antibody reads for single cell protein measurements.

# Important Tools:
After compilation tools are found in *./bin*

  - **PARSER**: Mapping of a barcode pattern to fastq-reads
  - **PROCESSING**: Generating a singleCell * ABCount matrix for those mapped reads
  
# Get started:

  `make install`  
  `make parser`  
  `make processing`  
 - Have a look in the Makefile for how to use tool
