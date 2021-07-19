#DEPENDANCIES: zlib, inoput is a ONE READ fastq file, there fore convert forward/ reverse fastqs into one e.g. with bbmap
install:
	#download and compile kseq
	mkdir Tools; cd ./Tools; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

parser:
	g++ src/Tools/AbFastqParser/FastqParser.cpp -o bin/parser -pthread -lz -lboost_program_options -I ./Tools/ --std=c++11

processing:
	g++ src/Tools/BarcodeProcessing/BarcodeProcessing.cpp -o bin/processing -pthread -lz -lboost_program_options -I ./Tools/ -I ./src/lib --std=c++11

run:
	./bin/parser -i ./src/test/test_data/inFastqTest.fastq -o ./output.csv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN][XXX][TTT] -m 1,2,1,1,1 -t 1 -b ./src/test/test_data/barcodeFile.txt
