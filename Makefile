#DEPENDANCIES: 
#	zlib /(input is a ONE READ fastq file, therefore convert forward/ reverse fastqs into one e.g. with fastq-join/)
#	boost

#STAR with a minimum version of  2.7.9, to also annotate gene ids (GX)
#WARNING:
# FOR NOW RNA-MAPPING IS NOT SUPPORTED UNDER WINDOWS: missing htslib installation!
#we need htslib to annotate mapped barcodes to the mapped genes in the BAM file

UNAME_S := $(shell uname -s)
VCPKG_ROOT ?= C:/vcpkg

CXXFLAGS = -g -Wall -Wextra	-Wsign-compare -O3

#system dependent boost flags
ifeq ($(UNAME_S),Linux)
    BOOST_FLAGS = -lboost_iostreams -lboost_program_options -lpthread -lz
	BOOST_INCLUDE =
	BOOST_LIB =
endif

ifeq ($(UNAME_S),Darwin)
    BOOST_FLAGS = -lboost_iostreams -lboost_program_options -lpthread -lz
	BOOST_INCLUDE =
	BOOST_LIB =
endif

# Check for any Windows-like environments: MINGW, MSYS, or CYGWIN
BOOST_LIB_NAMES = system thread program_options iostreams
ifneq (,$(findstring MINGW,$(UNAME_S))$(findstring MSYS,$(UNAME_S)))
    BOOST_INCLUDE = $(VCPKG_ROOT)/installed/x64-mingw-static/include
    BOOST_LIB = $(VCPKG_ROOT)/installed/x64-mingw-static/lib
	#dynamically detect the right boost suffix
	BOOST_SUFFIX := $(shell \
		for lib in $(BOOST_LIB)/libboost_system*.a; do \
			basename $$lib | sed -n 's/libboost_system\(.*\)\.a/\1/p'; \
		done | head -n 1)
    BOOST_FLAGS := -I$(BOOST_INCLUDE) -L$(BOOST_LIB) $(foreach lib,$(BOOST_LIB_NAMES),-lboost_$(lib)$(BOOST_SUFFIX)) -lpthread -lz -lwinpthread
endif

#only include boost flags if needed
BOOST_INCLUDE_FLAG := $(if $(BOOST_INCLUDE),-I$(BOOST_INCLUDE),)

#LDFLAGS
ifeq ($(OS), Windows_NT)
	LDFLAGS += -static -static-libgcc -static-libstdc++ -Wl,-Bstatic -lz -lwinpthread 
else ifeq ($(UNAME_S),Linux)
	LDFLAGS +=
else ifeq ($(UNAME_S),Darwin)
	LDFLAGS += 
endif

install:
	# download and compile kseq (we have a modified makefile to compile with rand on windows, which we move into the repo), 
	# we have a submodule edlib (git submodule add https://github.com/martinsos/edlib ./edlib;
	# no need to compile, we just add libraries and then compile with them), but we update it

	#create bin	
	mkdir bin

	#update submodules (edlib)
	git submodule update --init --recursive

	#we use seqtk to hadnel fastq files, we directly include the files and need to download them from git
	cd ./include; git clone https://github.com/lh3/seqtk --branch v1.3; mv Makefile ./seqtk/; mv rand_win.c ./seqtk/; cd ./seqtk; make;

	#install libboost for various systems LINUX/ WINDOWS/ macOS
	#TODO: we do not need all libboost-dev for LINUX and boost for macOS (check which libs are needed and install only those!)
	@if [ "$(UNAME_S)" = "Linux" ]; then \
		sudo apt-get update && sudo apt-get install -y zlib1g-dev libboost-all-dev libhts-dev; \
	elif echo "$(UNAME_S)" | grep -E -q "MINGW|MSYS"; then \
		vcpkg install zlib boost-asio boost-system boost-thread boost-iostreams boost-program-options --triplet x64-mingw-static; \
	elif [ "$(UNAME_S)" = "Darwin" ]; then \
		brew install zlib boost htslib; \
	fi

	#build htslib manually
	#cd ./include; git clone --recurse-submodules https://github.com/samtools/htslib.git; cd htslib; $(MAKE); $(MAKE) -C htslib install; cd ..

ESGI:
	make demultiplex
	make count
	

#parse fastq lines and map abrcodes to each sequence
demultiplex:
	g++ -c ./include/edlib/edlib/src/edlib.cpp -I ./include/edlib/edlib/include/ -I ./src/lib $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/lib/DemultiplexedStatistics.cpp -I ./include/ -I ./src/lib $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/lib/BarcodeMapping.cpp -I ./include/ -I ./src/lib $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/Demultiplexing/DemultiplexedResult.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/Demultiplexing/Demultiplexer.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/Demultiplexing/main.cpp -o main_demultiplex.o -I ./include/ -I ./src/lib -I src/tools/Demultiplexing $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ main_demultiplex.o DemultiplexedResult.o Demultiplexer.o BarcodeMapping.o DemultiplexedStatistics.o edlib.o -o ./bin/demultiplex $(LDFLAGS) $(BOOST_FLAGS)

#process the mapped sequences: correct for UMI-mismatches, then map barcodes to Protein, treatment, SinglecellIDs
count:
	g++ -c src/tools/FeatureCounting/BarcodeProcessingHandler.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/FeatureCounting/main.cpp -o main_count.o -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing $(BOOST_INCLUDE_FLAG) --std=c++17 $(CXXFLAGS)
	g++ main_count.o BarcodeProcessingHandler.o -o ./bin/count $(LDFLAGS) $(BOOST_FLAGS)

annotate:
	g++ -o ./bin/annotate src/tools/BarcodefileBamAnnotator/BarcodeBamAnnotator.cpp src/tools/BarcodefileBamAnnotator/main.cpp $(LDFLAGS) -lboost_iostreams -lboost_program_options -lhts

#a quality control tool: Mapping first Linker to whole sequence
demultiplexAroundLinker:
	g++ -c src/lib/BarcodeMapping.cpp -I ./include/ -I ./src/lib -I src/tools/Demultiplexing --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/DemultiplexAroundLinker/MappingAroundLinker.cpp -I ./include/ -I ./src/lib -I src/tools/DemultiplexAroundLinker --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/DemultiplexAroundLinker/main.cpp -I ./include/ -I ./src/tools/DemultiplexAroundLinker -I ./src/tools/FeatureCounting -I ./src/lib --std=c++17 $(CXXFLAGS)
	g++ main.o MappingAroundLinker.o BarcodeMapping.o -o ./bin/demultiplexAroundLinker -lpthread -lz -lboost_program_options -lboost_iostreams

testDemultiplexAroundLinker:
	#simple example with a few missing barcodes
	./bin/demultiplexAroundLinker -i ./src/test/test_data/inFile_DemAroundLinker.txt -o ./bin/output.tsv -p [CATGAGCGTCATG][NNNN][XXX][ATCAGTCAACAGATAAGCGA][NNNN] -m 1,1,1,1,1 -t 1 -b ./src/test/test_data/barcodeFile.txt
	(head -n 1 ./bin/DemultiplexedAroundLinker_output.tsv && tail -n +2 ./bin/DemultiplexedAroundLinker_output.tsv | LC_ALL=c sort)  > ./bin/demultiplexAroundLinker_output1_sorted.tsv
	diff ./src/test/test_data/demultiplexAroundLinker_output1_sorted.tsv ./bin/demultiplexAroundLinker_output1_sorted.tsv

	#more complex example with several barcodes missing in the middle that are all the same (sevel BCR missing with same barcodes)
	./bin/demultiplexAroundLinker -i ./src/test/test_data/inFile_DemAroundLinker_2.txt -o ./bin/output.tsv -p [CATGAGCGTCATG][NNNN][CATGAGCGTCATG][NNNN][XXX][ATCAGTCAACAGATAAGCGA][NNNN] -m 1,1,1,1,1,1,1 -t 1 -b ./src/test/test_data/barcodeFile_2.txt
	(head -n 1 ./bin/DemultiplexedAroundLinker_output.tsv && tail -n +2 ./bin/DemultiplexedAroundLinker_output.tsv | LC_ALL=c sort)  > ./bin/demultiplexAroundLinker_output2_sorted.tsv
	diff ./src/test/test_data/demultiplexAroundLinker_output2_sorted.tsv ./bin/demultiplexAroundLinker_output2_sorted.tsv

#Umiqual is a toll to analuze the quality of the CI reads based on the UMI. Imagine we have an explosion of barocode combinations
# we can use this tool to see in which BC round those combinations occure (based on the UMI)
umiqual:
	g++ -c src/tools/FeatureCounting/BarcodeProcessingHandler.cpp -I ./include/ -I ./src/lib -I ./src/tools/Demultiplexing --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/UmiQualityCheck/UmiQualityHelper.cpp -I ./include/ -I ./src/lib -I ./src/tools/FeatureCounting --std=c++17 $(CXXFLAGS)
	g++ -c src/tools/UmiQualityCheck/main.cpp -I ./include/ -I ./src/lib -I ./src/tools/FeatureCounting --std=c++17 $(CXXFLAGS)
	g++ main.o UmiQualityHelper.o BarcodeProcessingHandler.o -o ./bin/umiqual -lpthread -lz -lboost_program_options -lboost_iostreams

#this test is run on github.com
test:
	make demultiplex
	make test_demultiplex
	make bigTest
	make test_multipattern

	make count
	make test_count

test_detached:
	./bin/demultiplex -i ./src/test/test_data/test_detached/input_fw.fastq -r ./src/test/test_data/test_detached/input_rv.fastq -d 1 -o ./bin/ -p ./src/test/test_data/test_detached/patterns.txt -m ./src/test/test_data/test_detached/mismatches.txt -t 1 -n DETACHED -q 1 -f 1

test_multipattern:
	./bin/demultiplex -i ./src/test/test_data/test_multipatterns/input.txt -o ./bin/ -p ./src/test/test_data/test_multipatterns/patterns.txt -m ./src/test/test_data/test_multipatterns/mismatches.txt -t 1 -n MULTI -q 1 -f 1
	diff ./bin/MULTI_PATTERN1.tsv src/test/test_data/test_multipatterns/MULTI_PATTERN1.tsv
	diff ./bin/MULTI_PATTERN2.tsv src/test/test_data/test_multipatterns/MULTI_PATTERN2.tsv
	diff ./bin/MULTI_PATTERN3.tsv src/test/test_data/test_multipatterns/MULTI_PATTERN3.tsv

testUmiqual:
	./bin/umiqual -i ./src/test/test_data/testSet.txt.gz -o ./bin/processed_out.tsv -t 1 -b ./src/test/test_data/processingBarcodeFile.txt  -c 0,2,3,4 -a ./src/test/test_data/antibody.txt -x 1 -g ./src/test/test_data/treatment.txt -y 2 -u 2
	(head -n 1 ./bin/UmiQualityCheckprocessed_out.tsv && tail -n +2 ./bin/UmiQualityCheckprocessed_out.tsv | LC_ALL=c sort)  > ./bin/UmiQualityCheckprocessedSorted_out.tsv
	diff ./bin/UmiQualityCheckprocessedSorted_out.tsv ./src/test/test_data/UmiQualityCheckprocessed_out.tsv

#small test script for the parser, includes 1 perfect match, 6 matches with different types of mismatches below threshold, 2 mismatches above threshold
#and four mismatches due to barcodes that can not be uniquely identified
test_demultiplex_TODO:
	#test orde ron one thread
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -q true
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./bin/Demultiplexed_output.tsv
	diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./bin/StatsMismatches_output.tsv

	#test order with more threads
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 4 -b ./src/test/test_data/barcodeFile.txt
	(head -n 1 ./bin/Demultiplexed_output.tsv && tail -n +2 ./bin/Demultiplexed_output.tsv | LC_ALL=c sort)  > ./bin/DemultiplexedSorted_output.tsv
	(head -n 1 ./src/test/test_data/BarcodeMapping_output.tsv && tail -n +2 ./src/test/test_data/BarcodeMapping_output.tsv | LC_ALL=c sort)  > ./src/test/test_data/BarcodeMappingSorted_output.tsv
	diff ./src/test/test_data/BarcodeMappingSorted_output.tsv ./bin/DemultiplexedSorted_output.tsv
	#test paired end mapping
	./bin/demultiplexing -i ./src/test/test_data/smallTestPair_R1.fastq.gz -r ./src/test/test_data/smallTestPair_R2.fastq.gz -o ./bin/PairedEndTest -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 1,0,0,1,0,1,0,1,0,1 -t 1 -b ./src/test/test_data/processingBarcodeFile.txt
	diff ./bin/Demultiplexed_PairedEndTest ./src/test/test_data/result_pairedEnd
	#test paired end, where the overlapping barcode is a different one in forward/ reverse: example incldues also one read with forward barcode and reverse PhiX
	./bin/demultiplexing -i ./src/test/test_data/pairedtestR1.fastq -r ./src/test/test_data/pairedtestR2.fastq -o ./bin/Pairedtest2 -p [NNNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTCTCAAGCACGTGGAT][NNNNNNNN][AGTCGTACGCCGATGCGAAACATCGGCCAC][NNNNNNNN] -m 1,2,0,1,2,1,2,1,15,2 -t 1 -b ./src/test/test_data/processingBarcodefilewithStagger.txt
	diff ./bin/Demultiplexed_Pairedtest2 ./src/test/test_data/Demultiplexed_Pairedtest2.txt

	#test the offset for a next sequence after deletions at end of previous barcode
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest_4.fastq -o ./bin/output.tsv -p [NNNN][XXXX][ATATCAGTCGAAA][NNNN][AAAGCTCGATCAT] -m 1,1,2,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -q true
	(head -n 1 ./bin/Demultiplexed_output.tsv && tail -n +2 ./bin/Demultiplexed_output.tsv | LC_ALL=c sort)  > ./bin/DemultiplexedSorted_output.tsv
	diff ./src/test/test_data/DemultiplexedSortedElongationTest_output.tsv ./bin/DemultiplexedSorted_output.tsv

	#test the deletion of the barcode end (if a barcode is not mapped fully, the next one should have
	#variable starting points...)
	./bin/demultiplexing -i ./src/test/test_data/testDeletionsAtBarcodeEnd.fastq -o ./bin/testBarcodeDeletionsEnd.tsv -p [AGCTAGCTAAAA][TATA] -m 4,0 -t 1
	diff ./bin/Demultiplexed_testBarcodeDeletionsEnd.tsv ./src/test/test_data/result_barcodeDeletionEnd.tsv

	#test for double-UMI patterns
	./bin/demultiplexing -i ./src/test/test_data/inFastqDoubleUmiTest.fastq -o ./bin/testMultipleUmis.tsv -p [AAAA][XXXX][XXXX][TTTT] -m 1,1,1,1 -t 1
	diff ./bin/Demultiplexed_testMultipleUmis.tsv ./src/test/test_data/result_testMultipleUmis.tsv

	#paired end we do not find the UMIs with insertion/ deletion bcs from both sides we map different UMIs
	./bin/demultiplexing -i ./src/test/test_data/inFastqDoubleUmiTest_1.fastq -r ./src/test/test_data/inFastqDoubleUmiTest_2.fastq -o ./bin/testMultipleUmis_PairedEnd.tsv -p [AAAA][XXXX][XXXX][TTTT] -m 1,1,1,1 -t 1
	diff ./bin/Demultiplexed_testMultipleUmis_PairedEnd.tsv ./src/test/test_data/result_testMultipleUmis_PairedEnd.tsv

test_demultiplex:
	#test order on one thread
	./bin/demultiplex -i ./src/test/test_data/inFastqTest.fastq -o ./bin/ -p ./src/test/test_data/test1Pattern.txt -m ./src/test/test_data/test1MM.txt -t 1 -n TEST -q 1
	diff ./src/test/test_data/BarcodeMapping_output.tsv ./bin/TEST_TEST1.tsv
	
	#statistics not imlemented yet for new alignment method ...
	#diff ./src/test/test_data/StatsBarcodeMappingErrors_output.tsv ./bin/StatsMismatches_output.tsv

#sometimes several barcodes can encode for the same cell (e.g., look at SIGNALseq where two different barcodes tag
#poly-A and randomHexamer reads with two different barcodes), we can tell the 'count' tool to collapse those SC-barcodes
test_barcode_merging:
	./bin/count -i ./src/test/test_data/test_barcodeMerging/demultiplexedReads.tsv -o ./bin/barcodeMergingCounts.tsv -t 2 -d ./src/test/test_data/test_barcodeMerging -c 0,1 -a ./src/test/test_data/test_barcodeMerging/antibody.txt -x 3 -u 2 -m 1 -w ./src/test/test_data/test_barcodeMerging/mergeBarcodes.tsv -s 1
	(head -n 1 ./bin/ABbarcodeMergingCounts.tsv && tail -n +2 ./bin/ABbarcodeMergingCounts.tsv | LC_ALL=c sort) > ./bin/sortedABbarcodeMergingCounts.tsv
	diff ./src/test/test_data/test_barcodeMerging/ABbarcodeMergingCounts.tsv ./bin/sortedABbarcodeMergingCounts.tsv

#TO DO:
move_this_to_test_ezgi:
	#test order with more threads
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest.fastq -o ./bin/output.tsv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,1,2 -t 4 -b ./src/test/test_data/barcodeFile.txt
	(head -n 1 ./bin/Demultiplexed_output.tsv && tail -n +2 ./bin/Demultiplexed_output.tsv | LC_ALL=c sort)  > ./bin/DemultiplexedSorted_output.tsv
	(head -n 1 ./src/test/test_data/BarcodeMapping_output.tsv && tail -n +2 ./src/test/test_data/BarcodeMapping_output.tsv | LC_ALL=c sort)  > ./src/test/test_data/BarcodeMappingSorted_output.tsv
	diff ./src/test/test_data/BarcodeMappingSorted_output.tsv ./bin/DemultiplexedSorted_output.tsv
	#test paired end mapping
	./bin/demultiplexing -i ./src/test/test_data/smallTestPair_R1.fastq.gz -r ./src/test/test_data/smallTestPair_R2.fastq.gz -o ./bin/PairedEndTest -p [NNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTACGC][NNNNNNNN][CGAAGTCGTACGCCGATG][NNNNNNNN] -m 1,0,0,1,0,1,0,1,0,1 -t 1 -b ./src/test/test_data/processingBarcodeFile.txt
	diff ./bin/Demultiplexed_PairedEndTest ./src/test/test_data/result_pairedEnd
	#test paired end, where the overlapping barcode is a different one in forward/ reverse: example incldues also one read with forward barcode and reverse PhiX
	./bin/demultiplexing -i ./src/test/test_data/pairedtestR1.fastq -r ./src/test/test_data/pairedtestR2.fastq -o ./bin/Pairedtest2 -p [NNNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTCTCAAGCACGTGGAT][NNNNNNNN][AGTCGTACGCCGATGCGAAACATCGGCCAC][NNNNNNNN] -m 1,2,0,1,2,1,2,1,15,2 -t 1 -b ./src/test/test_data/processingBarcodefilewithStagger.txt
	diff ./bin/Demultiplexed_Pairedtest2 ./src/test/test_data/Demultiplexed_Pairedtest2.txt

	#test the offset for a next sequence after deletions at end of previous barcode
	./bin/demultiplexing -i ./src/test/test_data/inFastqTest_4.fastq -o ./bin/output.tsv -p [NNNN][XXXX][ATATCAGTCGAAA][NNNN][AAAGCTCGATCAT] -m 1,1,2,1,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -q true
	(head -n 1 ./bin/Demultiplexed_output.tsv && tail -n +2 ./bin/Demultiplexed_output.tsv | LC_ALL=c sort)  > ./bin/DemultiplexedSorted_output.tsv
	diff ./src/test/test_data/DemultiplexedSortedElongationTest_output.tsv ./bin/DemultiplexedSorted_output.tsv

	#test the deletion of the barcode end (if a barcode is not mapped fully, the next one should have
	#variable starting points...)
	./bin/demultiplexing -i ./src/test/test_data/testDeletionsAtBarcodeEnd.fastq -o ./bin/testBarcodeDeletionsEnd.tsv -p [AGCTAGCTAAAA][TATA] -m 4,0 -t 1
	diff ./bin/Demultiplexed_testBarcodeDeletionsEnd.tsv ./src/test/test_data/result_barcodeDeletionEnd.tsv

	#test for double-UMI patterns
	./bin/demultiplexing -i ./src/test/test_data/inFastqDoubleUmiTest.fastq -o ./bin/testMultipleUmis.tsv -p [AAAA][XXXX][XXXX][TTTT] -m 1,1,1,1 -t 1
	diff ./bin/Demultiplexed_testMultipleUmis.tsv ./src/test/test_data/result_testMultipleUmis.tsv

	#paired end we do not find the UMIs with insertion/ deletion bcs from both sides we map different UMIs
	./bin/demultiplexing -i ./src/test/test_data/inFastqDoubleUmiTest_1.fastq -r ./src/test/test_data/inFastqDoubleUmiTest_2.fastq -o ./bin/testMultipleUmis_PairedEnd.tsv -p [AAAA][XXXX][XXXX][TTTT] -m 1,1,1,1 -t 1
	diff ./bin/Demultiplexed_testMultipleUmis_PairedEnd.tsv ./src/test/test_data/result_testMultipleUmis_PairedEnd.tsv

test_count:
#origional first test with several basic examples
	./bin/count -i ./src/test/test_data/testSet.txt -o ./bin/processed_out.tsv -t 2 -d ./src/test/test_data -c 0,5,7,9 -a ./src/test/test_data/antibody.txt -x 3 -g ./src/test/test_data/treatment.txt -y 5 -u 2 -f 0.9
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	diff ./src/test/test_data/sortedABprocessed_out.tsv ./bin/sortedABprocessed_out.tsv

#test with multiple UMIs
	./bin/count -i ./src/test/test_data/testTwoUMIs.txt -o ./bin/2UMIs_out.tsv -t 2 -d ./src/test/test_data -c 0,5,7,9 -a ./src/test/test_data/antibody.txt -x 3 -g ./src/test/test_data/treatment.txt -y 5 -u 2,10 -f 0.9 -z 0
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	diff ./src/test/test_data/sortedABprocessedTwoUMIs_out.tsv ./bin/sortedABprocessed_out.tsv

#test processing of the barcodes, includes several UMIs with mismatches, test the mapping of barcodes to unique CellIDs, ABids, treatments
testProcessing:
#testing the removal of one wrong read bcs of different AB-Sc for same UMI
	./bin/count -i ./src/test/test_data/testSet_2.txt.gz -o ./bin/processed_out.tsv -t 2 -d ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,5 -a ./src/test/test_data/antibody_2.txt -x 3 -g ./src/test/test_data/treatment_2.txt -y 5 -u 2 -f 0.9
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedABprocessed_out.tsv ./src/test/test_data/sortedABprocessed_2_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/sortedUMIprocessed_2_out.tsv
#testing removal of two reads bcs both have different treatments for same SC
	./bin/processing -i ./src/test/test_data/test_treatmentReadRemoval.txt.gz -o ./bin/processed_out.tsv -t 2 -d ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,2 -a ./src/test/test_data/antibody_2.txt -x 1 -g ./src/test/test_data/treatment_2.txt -y 2 -u 2 -f 0.9
	(head -n 1 ./bin/ABprocessed_out.tsv && tail -n +2 ./bin/ABprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedABprocessed_out.tsv
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedABprocessed_out.tsv ./src/test/test_data/sortedABprocessed_treatment_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/sortedUMIprocessed_treatment_out.tsv
#testing the specific removal of two lines bcs the same SC (with two lines only) has two different treatments (additional one line is removed for UMI)
	./bin/processing -i ./src/test/test_data/test_treatmentReadRemoval_Log.txt.gz -o ./bin/processed_out.tsv -t 2 -d ./src/test/test_data/processingBarcodeFile_2.txt  -c 0 -a ./src/test/test_data/antibody_2.txt -x 1 -g ./src/test/test_data/treatment_2.txt -y 2 -u 2 -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt -f 0.9
	(head -n 1 ./bin/LOGprocessed_out.tsv && tail -n +2 ./bin/LOGprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedLOGprocessed_out.tsv
	diff ./bin/sortedLOGprocessed_out.tsv ./src/test/test_data/sortedLOGprocessed_treatment_out.tsv
#test EditDist for UMIs
	./bin/processing -i ./src/test/test_data/umiEditDistTest.txt.gz -o ./bin/processed_out.tsv -t 2 -d ./src/test/test_data/processingBarcodeFile_2.txt  -c 0,2 -a ./src/test/test_data/antibody_2.txt -x 1 -g ./src/test/test_data/treatment_2.txt -y 2 -u 2 -f 0.9
	(head -n 1 ./bin/UMIprocessed_out.tsv && tail -n +2 ./bin/UMIprocessed_out.tsv | LC_ALL=c sort) > ./bin/sortedUMIprocessed_out.tsv
	diff ./bin/sortedUMIprocessed_out.tsv ./src/test/test_data/UMIprocessed_out_editTest.tsv

#testing the whole analysis pipeline to smoothly run through with a few additional test scenarios
testAnalysis:
#a basic test from mostly already existing files, just to check tool runs through
	php ./src/Pipelines/analyze.php -a ./src/test/test_data/antibody_3.txt -i ./src/test/test_data/inFastqTest_2.fastq -o ./bin/AnalysisTestOutput -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][GATCAT] -m 1,4,1,2,2 -t 1 -b ./src/test/test_data/barcodeFile.txt -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt -x 1 -c 0 -f 0.9 -h true
	(head -n 1 ./bin/AnalysisTestOutput/ABProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/ABProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv ./src/test/test_data/ABProcessing_AnalyzePipeline.tsv
#more elaborate test with common real world scenarios: a test with three single cells with two different proteins each, one cell beeing removed bcs it has no guide read, two reads removed for same UMI but different proteins, and one cell
#beeing removed for having different guide reads, also one duplicate UMI read is present that has to only be counted ONCE
	php ./src/Pipelines/analyze.php -a ./src/test/test_data/antibody_3.txt -i ./src/test/test_data/inFastqTest_3.fastq -o ./bin/AnalysisTestOutput -p [NNNN][ATCAGTCA][NNNN][ACAGATAAGCGA][NNNN][XXXX] -m 1,2,1,2,1,1 -t 1 -b ./src/test/test_data/barcodeFile_2.txt -g ./src/test/test_data/guideTest_class_seqs.txt -n ./src/test/test_data/guideTest_class_names.txt -x 2 -c 0,1 -f 0.9 -h true
	(head -n 1 ./bin/AnalysisTestOutput/ABProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/ABProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/ABProcessing_Sorted.tsv ./src/test/test_data/ABProcessing_AnalyzePipeline_2.tsv
#also test the output of UMI counts (so was not tested anywhere else, so just added it here out of convenience) tests umi duplicates without mismatches and
#also one umi that has one deletion but should still be counted twice as we allow for one MM
	(head -n 1 ./bin/AnalysisTestOutput/UMIProcessing.tsv && tail -n +2 ./bin/AnalysisTestOutput/UMIProcessing.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/UMIProcessing_Sorted.tsv
	diff ./bin/AnalysisTestOutput/UMIProcessing_Sorted.tsv ./src/test/test_data/UmiProcessed_Test.tsv
#test combination of running demultiplexing/ and then processing on actual data (duplicated some reads to make sure they r not counted twice and dublicated two reads and added a new UMI to make sure they r counted)
	./bin/demultiplexing -i ./src/test/test_data/testFullAnalysisR1.fastq -r ./src/test/test_data/testFullAnalysisR2.fastq -o ./bin/AnalysisTestOutput/FullAnalysis.tsv -p [NNNNNNNNN][CTTGTGGAAAGGACGAAACACCG][XXXXXXXXXXXXXXX][NNNNNNNNNN][GTTTTAGAGCTAGAAATAGCAA][NNNNNNNN][CGAATGCTCTGGCCTCTCAAGCACGTGGAT][NNNNNNNN][AGTCGTACGCCGATGCGAAACATCGGCCAC][NNNNNNNN] -b ./src/test/test_data/barcodesFullAnalysis.txt -m 1,15,0,1,15,1,15,1,15,1 -t 1 -q true -f true
	gzip ./bin/AnalysisTestOutput/Demultiplexed_FullAnalysis.tsv
	./bin/processing -i ./bin/AnalysisTestOutput/Demultiplexed_FullAnalysis.tsv.gz -o ./bin/AnalysisTestOutput/FullAnalysis_ABCOUNT_RESULT.tsv -b ./src/test/test_data/barcodesFullAnalysis.txt -a ./src/test/test_data/antibodiesFullAnalysis.txt -x 1 -c 0,2,3,4 -u 0 -t 1 -d ./src/test/test_data/treatmentsFullAnalysis.txt -y 2
	rm ./bin/AnalysisTestOutput/Demultiplexed_FullAnalysis.tsv.gz
	(head -n 1 ./bin/AnalysisTestOutput/ABFullAnalysis_ABCOUNT_RESULT.tsv && tail -n +2 ./bin/AnalysisTestOutput/ABFullAnalysis_ABCOUNT_RESULT.tsv | LC_ALL=c sort) > ./bin/AnalysisTestOutput/ABFullAnalysis_ABCOUNT_RESULT_SORTED.tsv
	diff ./bin/AnalysisTestOutput/ABFullAnalysis_ABCOUNT_RESULT_SORTED.tsv ./src/test/test_data/FullAnalysis_ABCOUNT_RESULT.tsv


debug:
	time ./bin/demultiplex -i ./test1.fastq.gz -r ./test2.fastq.gz -o ./bin/ -n debug -p ./CITestData/background_data/pattern.txt -m ./CITestData/background_data/mismatches.txt -t 1 -f 1 -q 1


#some of Kathys CI data with GUIDE and PROTEIN
CITest:
	./bin/demultiplex -i ./CITestData/CITest_1.fastq.gz -r ./CITestData/CITest_2.fastq.gz -o ./bin/ -n CITEST -p ./CITestData/background_data/pattern.txt -m ./CITestData/background_data/mismatches.txt -t 10 -f 1 -q 1

#single AB pattern
bigTest:
	#time ./bin/demultiplex -i tmp.fastq -o ./bin/ -p ./src/test/test_data/test_input/barcodePatternsBig.txt -m ./src/test/test_data/test_input/barcodeMismatchesBig.txt -t 1 -f 1 -q 1
	./bin/demultiplex -i ./src/test/test_data/test_input/testBig.fastq.gz -o ./bin/ -p ./src/test/test_data/test_input/barcodePatternsBig.txt -m ./src/test/test_data/test_input/barcodeMismatchesBig.txt -t 1 -f 1 -q 1

bigTest2:
	time ./bin/demultiplex -i ./src/test/test_data/test_input/testBig2.fastq.gz -o ./bin/ -p ./src/test/test_data/test_input/barcodePatternsBig.txt -m /DATA/t.stohn/SCDemultiplexing/src/test/test_data/test_input/barcodeMismatchesBig.txt -t 50 -f 1
bigTest3:
	time ./bin/demultiplex -i ./src/test/test_data/test_input/testBig3.fastq.gz -o ./bin/ -p /DATA/t.stohn/SCDemultiplexing/src/test/test_data/test_input/barcodePatternsBig.txt -m /DATA/t.stohn/SCDemultiplexing/src/test/test_data/test_input/barcodeMismatchesBig.txt -t 70 -f 1

#makes no sense since we have only forward reads...
#make a small test and use fw and rv files
bigTestRNA:
	time ./bin/demultiplex -i ./src/test/test_data/test_input/testBig.fastq.gz -o ./bin/ -p ./src/test/test_data/test_input/barcodePatternsBigRNA.txt -m ./src/test/test_data/test_input/barcodeMismatchesBigRNA.txt -t 10 -f 1

testGuideMapping:
	./bin/demultiplex -i ./src/test/test_data/guideTestInput.txt -o ./bin/output.tsv -p [NNNNNNN][GTTTAAA][XXXXXXXXXX][NNNNNNNNNN] -m 1,1,1,1 -t 1 -b ./src/test/test_data/barcodesGuideTest.txt -c ./src/test/test_data/guidesGuideTest.txt -e 1
	diff ./bin/Demultiplexed_guideReadsoutput.tsv ./src/test/test_data/guideTestGUIDEOutput.txt
	diff ./bin/Demultiplexed_output.tsv ./src/test/test_data/guideTestABOutput.txt
