#DEPENDANCIES: zlib, input is a ONE READ fastq file, therefore convert forward/ reverse fastqs into one e.g. with fastq-join
install:
	#download and compile kseq
	mkdir include; cd ./include; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

#parse fastq lines and map abrcodes to each sequence
demultiplexing:
	g++ -c src/lib/BarcodeMapping.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/Demultiplexing/DemultiplexedLinesWriter.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/Demultiplexing/main.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing --std=c++17
	g++ main.o DemultiplexedLinesWriter.o BarcodeMapping.o -o ./bin/demultiplexing -lpthread -lz -lboost_program_options -lboost_iostreams

#a quality control tool: Mapping first Linker to whole sequence
demultiplexAroundLinker:
	g++ -c src/lib/BarcodeMapping.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/BarcodeProcessing/UmiDataParser.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/DemultiplexAroundLinker/main.cpp -I ./include/ -I ./src/tools/Demultiplexing -I ./src/tools/BarcodeProcessing -I ./src/lib --std=c++17
	g++ main.o UmiDataParser.o BarcodeMapping.o -o ./bin/demultiplexAroundLinker -lpthread -lz -lboost_program_options -lboost_iostreams

#process the mapped sequences: correct for UMI-mismatches, then map barcodes to Protein, treatment, SinglecellIDs
processing:
	g++ -c src/tools/BarcodeProcessing/UmiDataParser.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/BarcodeProcessing/main.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ main.o UmiDataParser.o -o ./bin/processing -lpthread -lz -lboost_program_options -lboost_iostreams

#small test script for the parser, includes 1 perfect match, 6 matches with different types of mismatches below threshold, 2 mismatches above threshold
#and four mismatches due to barcodes that can not be uniquely identified
testDemultiplexing:
	#test orde ron one thread
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -q true
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./bin/BarcodeMapping_output.tsv
	diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./bin/StatsBarcodeMappingErrors_output.tsv

	#test order with more threads
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 4 -b ./src/test/test_data/barcodeFile.txt
	(head -n 1 ./bin/BarcodeMapping_output.tsv && tail -n +2 ./bin/BarcodeMapping_output.tsv | sort)  > ./bin/BarcodeMappingSorted_output.tsv
	(head -n 1 ./src/test/test_data/BarcodeMapping_output.tsv && tail -n +2 ./src/test/test_data/BarcodeMapping_output.tsv | sort)  > ./src/test/test_data/BarcodeMappingSorted_output.tsv

	diff ./src/test/test_data/BarcodeMappingSorted_output.tsv ./bin/BarcodeMappingSorted_output.tsv

#test processing of the barcodes, includes several UMIs with mismatches, test the mapping of barcodes to unique CellIDs, ABids, treatments
testProcessing:
	./bin/processing -i ./src/test/test_data/testSet.txt.gz -o ./bin/processed_out.tsv -t 1 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -g ./src/test/test_data/treatment.txt -y 2 -m 2
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./src/test/test_data/UMIprocessed_out ./bin/sortedUMIprocessed_out.tsv 
	diff ./src/test/test_data/STATSprocessed_out ./bin/STATSprocessed_out.tsv
	diff ./src/test/test_data/ABprocessed_out ./bin/ABprocessed_out.tsv

testAnalysis:
	php ./src/Pipelines/analyze.php -i ./src/test/test_data/inFastqTest.fastq -o ./bin/AnalysisTestOutput.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 1

bigTest:
	./bin/demultiplexing -i ./src/test/test_data/test2000fastq.gz -o ./bin/output.tsv -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 7,13,0,8,13,6,13,4,13,4 -t 5 -b ./src/test/test_data/processingBarcodeFile.txt

bigFailure:
	/Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/bin/parserQC -t 1 -f /Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/bin/FailedLines_output.tsv  -o FAILURE.txt -b /Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/src/test/test_data/processingBarcodeFile.txt -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 1,3,0,1,3,1,3,1,3,1