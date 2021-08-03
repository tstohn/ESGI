#DEPENDANCIES: zlib, input is a ONE READ fastq file, therefore convert forward/ reverse fastqs into one e.g. with fastq-join
install:
	#download and compile kseq
	mkdir Tools; cd ./Tools; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

#parse fastq lines and map abrcodes to each sequence
parser:
	g++ src/Tools/AbFastqParser/FastqParser.cpp -o bin/parser -pthread -lz -lboost_program_options -I ./Tools/ --std=c++17

#process the mapped sequences: correct for UMI-mismatches, then map barcodes to Protein, treatment, SinglecellIDs
processing:
	g++ -c src/Tools/BarcodeProcessing/UmiDataParser.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ -c src/Tools/BarcodeProcessing/BarcodeProcessing.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ BarcodeProcessing.o UmiDataParser.o -o ./bin/processing -lpthread -lz -lboost_program_options -lboost_iostreams

#small test script for the parser, includes 1 perfect match, 6 matches with different types of mismatches below threshold, 2 mismatches above threshold
#and four mismatches due to barcodes that can not be uniquely identified
testParser:
	./bin/parser -i ./src/test/test_data/inFastqTest.fastq -o ./Output/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][TTT] -m 1,2,1,1,1 -t 1 -b ./src/test/test_data/barcodeFile.txt
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./Output/BarcodeMapping_output.tsv
	diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./Output/StatsBarcodeMappingErrors_output.tsv

#test processing of the barcodes, includes several UMIs with mismatches, test the mapping of barcodes to unique CellIDs, ABids, treatments
testProcessing:
	./bin/processing -i ./src/test/test_data/testSet.txt.gz -o ./Output/processed_out.tsv -t 1 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -g ./src/test/test_data/treatment.txt -y 2 -m 2
	diff ./src/test/test_data/UMIprocessed_out ./Output/UMIprocessed_out.tsv 
	diff ./src/test/test_data/STATSprocessed_out ./Output/STATSprocessed_out.tsv
	diff ./src/test/test_data/ABprocessed_out ./Output/ABprocessed_out.tsv
