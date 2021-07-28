#DEPENDANCIES: zlib, inoput is a ONE READ fastq file, there fore convert forward/ reverse fastqs into one e.g. with bbmap
install:
	#download and compile kseq
	mkdir Tools; cd ./Tools; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

parser:
	g++ src/Tools/AbFastqParser/FastqParser.cpp -o bin/parser -pthread -lz -lboost_program_options -I ./Tools/ --std=c++17

processing:
	g++ -c src/Tools/BarcodeProcessing/UmiDataParser.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ -c src/Tools/BarcodeProcessing/BarcodeProcessing.cpp -I ./Tools/ -I ./src/lib -I ./src/Tools/AbFastqParser --std=c++17
	g++ BarcodeProcessing.o UmiDataParser.o -o ./bin/processing_fast -lpthread -lz -lboost_program_options -lboost_iostreams

run:
	./bin/parser -i ./src/test/test_data/inFastqTest.fastq -o ./output.csv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][TTT] -m 1,2,1,1,1 -t 1 -b ./src/test/test_data/barcodeFile.txt

run_processing:
	./bin/processing_fast -i /Users/t.stohn/Desktop/KATHY/TestDataRun/barcodeTest3First1000.csv.gz -o ../../../KATHY/CIAB_OUTPUT/processed_out -t 1 -b /Users/t.stohn/Desktop/KATHY/TestDataRun/barcodeFile.txt  -c 0,2,3,4 -a /Users/t.stohn/Desktop/KATHY/TestDataRun/antibody.txt -x 1 -g /Users/t.stohn/Desktop/KATHY/TestDataRun/treatment.txt -y 2 -m 4
