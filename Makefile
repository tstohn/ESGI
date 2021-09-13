#DEPENDANCIES: zlib, input is a ONE READ fastq file, therefore convert forward/ reverse fastqs into one e.g. with fastq-join
install:
	#download and compile kseq
	mkdir Tools; cd ./Tools; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

#parse fastq lines and map abrcodes to each sequence
parser:
	g++ src/Tools/AbFastqParser/FastqParser.cpp -o bin/parser -pthread -lz -lboost_program_options -I ./Tools/ -I ./src/lib --std=c++17

qualityControl:
	g++ -c src/Tools/BarcodeProcessing/UmiDataParser.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ -c src/Tools/AbFastqParserQC/FastqParserQC.cpp -I ./Tools/ -I ./src/Tools/AbFastqParser -I ./src/Tools/BarcodeProcessing -I ./src/lib --std=c++17
	g++ FastqParserQC.o UmiDataParser.o -o ./bin/parserQC -lpthread -lz -lboost_program_options -lboost_iostreams

#process the mapped sequences: correct for UMI-mismatches, then map barcodes to Protein, treatment, SinglecellIDs
processing:
	g++ -c src/Tools/BarcodeProcessing/UmiDataParser.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ -c src/Tools/BarcodeProcessing/BarcodeProcessing.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ BarcodeProcessing.o UmiDataParser.o -o ./bin/processing -lpthread -lz -lboost_program_options -lboost_iostreams

#small test script for the parser, includes 1 perfect match, 6 matches with different types of mismatches below threshold, 2 mismatches above threshold
#and four mismatches due to barcodes that can not be uniquely identified
testParser:
	./bin/parser -i ./src/test/test_data/inFastqTest.fastq -o ./Output/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./Output/BarcodeMapping_output.tsv
	diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./Output/StatsBarcodeMappingErrors_output.tsv

#test processing of the barcodes, includes several UMIs with mismatches, test the mapping of barcodes to unique CellIDs, ABids, treatments
testProcessing:
	./bin/processing -i ./src/test/test_data/testSet.txt.gz -o ./Output/processed_out.tsv -t 1 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -g ./src/test/test_data/treatment.txt -y 2 -m 2
	(head -n 1 ./Output/UMIprocessed_out.tsv && tail -n +2 ./Output/UMIprocessed_out.tsv | sort) > ./Output/sortedUMIprocessed_out.tsv
	diff ./src/test/test_data/UMIprocessed_out ./Output/sortedUMIprocessed_out.tsv 
	diff ./src/test/test_data/STATSprocessed_out ./Output/STATSprocessed_out.tsv
	diff ./src/test/test_data/ABprocessed_out ./Output/ABprocessed_out.tsv

bigTest:
	./bin/parser -i ./src/test/test_data/test2000fastq.gz -o ./Output/output.tsv -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 4,10,4,4,10,4,10,4,10,4 -t 5 -b /Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/src/test/test_data/processingBarcodeFile.txt

