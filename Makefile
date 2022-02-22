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
	g++ -c src/tools/BarcodeProcessing/BarcodeProcessingHandler.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/BarcodeProcessing/main.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ main.o BarcodeProcessingHandler.o -o ./bin/processing -lpthread -lz -lboost_program_options -lboost_iostreams

umiqual:
	g++ -c src/tools/BarcodeProcessing/BarcodeProcessingHandler.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17
	g++ -c src/tools/UmiQualityCheck/UmiQualityHelper.cpp -I ./include/ -I ./src/lib -I ./src/tools/BarcodeProcessing --std=c++17
	g++ -c src/tools/UmiQualityCheck/main.cpp -I ./include/ -I ./src/lib -I ./src/tools/BarcodeProcessing --std=c++17
	g++ main.o UmiQualityHelper.o BarcodeProcessingHandler.o -o ./bin/umiqual -lpthread -lz -lboost_program_options -lboost_iostreams

testUmiqual:
	./bin/umiqual -i ./src/test/test_data/testSet.txt.gz -o ./bin/processed_out.tsv -t 1 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -g ./src/test/test_data/treatment.txt -y 2 -u 2
	diff ./bin/UmiQualityCheckprocessed_out.tsv ./src/test/test_data/UmiQualityCheckprocessed_out.tsv

#small test script for the parser, includes 1 perfect match, 6 matches with different types of mismatches below threshold, 2 mismatches above threshold
#and four mismatches due to barcodes that can not be uniquely identified
testDemultiplexing:
	#test orde ron one thread
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -q true
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./bin/Demultiplexed_output.tsv
	diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./bin/StatsMismatches_output.tsv

	#test order with more threads
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 4 -b ./src/test/test_data/barcodeFile.txt
	(head -n 1 ./bin/Demultiplexed_output.tsv && tail -n +2 ./bin/Demultiplexed_output.tsv | sort)  > ./bin/DemultiplexedSorted_output.tsv
	(head -n 1 ./src/test/test_data/BarcodeMapping_output.tsv && tail -n +2 ./src/test/test_data/BarcodeMapping_output.tsv | sort)  > ./src/test/test_data/BarcodeMappingSorted_output.tsv
	diff ./src/test/test_data/BarcodeMappingSorted_output.tsv ./bin/DemultiplexedSorted_output.tsv
	#test paired end mapping
	./bin/demultiplexing -i ./src/test/test_data/smallTestPair_R1.fastq.gz -r ./src/test/test_data/smallTestPair_R2.fastq.gz -o ./bin/PairedEndTest -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 1,0,0,1,0,1,0,1,0,1 -t 1 -b ./src/test/test_data/processingBarcodeFile.txt
	diff ./bin/Demultiplexed_PairedEndTest ./src/test/test_data/result_pairedEnd

#test processing of the barcodes, includes several UMIs with mismatches, test the mapping of barcodes to unique CellIDs, ABids, treatments
testProcessing:
#origional first test with several basic examples
	./bin/processing -i ./src/test/test_data/testSet.txt.gz -o ./bin/processed_out.tsv -t 2 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -d ./src/test/test_data/treatment.txt -y 2 -u 2
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	diff ./src/test/test_data/sortedABprocessed_out.tsv ./bin/sortedABprocessed_out.tsv
#testing the removal of one wrong read bcs of different AB-Sc for same UMI
	./bin/processing -i ./src/test/test_data/testSet_2.txt.gz -o ./bin/processed_out.tsv -t 2 -b ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,2 -a ./src/test/test_data/antibody_2.txt -x 1 -d ./src/test/test_data/treatment_2.txt -y 2 -u 2
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedABprocessed_out.tsv ./src/test/test_data/sortedABprocessed_2_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/sortedUMIprocessed_2_out.tsv
#testing removal of two reads bcs both have different treatments for same SC
	./bin/processing -i ./src/test/test_data/test_treatmentReadRemoval.txt.gz -o ./bin/processed_out.tsv -t 2 -b ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,2 -a ./src/test/test_data/antibody_2.txt -x 1 -d ./src/test/test_data/treatment_2.txt -y 2 -u 2
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedABprocessed_out.tsv ./src/test/test_data/sortedABprocessed_treatment_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/sortedUMIprocessed_treatment_out.tsv
#testing the specific removal of two lines bcs the same SC (with two lines only) has two different treatments (additional one line is removed for UMI)
	./bin/processing -i ./src/test/test_data/test_treatmentReadRemoval_Log.txt.gz -o ./bin/processed_out.tsv -t 2 -b ./src/test/test_data/processingBarcodeFile_2.txt  -c 0 -a ./src/test/test_data/antibody_2.txt -x 1 -d ./src/test/test_data/treatment_2.txt -y 2 -u 2 -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt
	(head -n 1 ./bin/LOGprocessed_out.tsv && tail -n +2 ./bin/LOGprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedLOGprocessed_out.tsv
	diff ./bin/sortedLOGprocessed_out.tsv ./src/test/test_data/sortedLOGprocessed_treatment_out.tsv
#test EditDist for UMIs
	./bin/processing -i ./src/test/test_data/umiEditDistTest.txt.gz -o ./bin/processed_out.tsv -t 2 -b ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,2 -a ./src/test/test_data/antibody_2.txt -x 1 -d ./src/test/test_data/treatment_2.txt -y 2 -u 2
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/UMIprocessed_out_editTest.tsv

#testing the whole analysis pipeline to smoothly run through with a few additional test scenarios
testAnalysis:
#a basic test from mostly already existing files, just to check tool runs through
	php ./src/Pipelines/analyze.php -a ./src/test/test_data/antibody_3.txt -i ./src/test/test_data/inFastqTest_2.fastq -o ./bin/AnalysisTestOutput -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,0,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt -x 1 -c 0 
	(head -n 1 ./bin/AnalysisTestOutput/ABProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/ABProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv ./src/test/test_data/ABProcessing_AnalyzePipeline.tsv
#more elaborate test with common real world scenarios: a test with three single cells with two different proteins each, one cell beeing removed bcs it has no guide read, two reads removed for same UMI but different proteins, and one cell
#beeing removed for having different guide reads, also one duplicate UMI read is present that has to only be counted ONCE
	php ./src/Pipelines/analyze.php -a ./src/test/test_data/antibody_3.txt -i ./src/test/test_data/inFastqTest_3.fastq -o ./bin/AnalysisTestOutput -p [NNNN][ATCAGTCA][NNNN][ACAGATAAGCGA][NNNN][XXXX] -m 1,2,1,2,1,1 -t 1 -b ./src/test/test_data/barcodeFile_2.txt -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt -x 2 -c 0,1
	(head -n 1 ./bin/AnalysisTestOutput/ABProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/ABProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv ./src/test/test_data/ABProcessing_AnalyzePipeline_2.tsv
#also test the output of UMI counts (so was not tested anywhere else, so just added it here out of convenience) tests umi duplicates without mismatches and
#also one umi that has one deletion but should still be counted twice as we allow for one MM
	(head -n 1 ./bin/AnalysisTestOutput/UMIProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/UMIProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/UMIProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/UMIProcessing_Sorted.tsv ./src/test/test_data/UmiProcessed_Test.tsv


bigTest:
	./bin/demultiplexing -i ./src/test/test_data/test2000fastq.gz -o ./bin/output.tsv -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 7,13,0,8,13,6,13,4,13,4 -t 5 -b ./src/test/test_data/processingBarcodeFile.txt

bigFailure:
	/Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/bin/parserQC -t 1 -f /Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/bin/FailedLines_output.tsv  -o FAILURE.txt -b /Users/t.stohn/Desktop/Normalization/PIPELINE/SingleCellGenomeTools/src/test/test_data/processingBarcodeFile.txt -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 1,3,0,1,3,1,3,1,3,1