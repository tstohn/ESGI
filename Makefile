#DEPENDANCIES: zlib, inoput is a ONE READ fastq file, there fore convert forward/ reverse fastqs into one e.g. with bbmap
install:
	#download and compile kseq
	mkdir Tools; cd ./Tools; git clone https://github.com/lh3/seqtk --branch v1.3; cd ./seqtk; make
	mkdir bin

parser:
	g++ src/Tools/FastqParser.cpp -o bin/parser -pthread -lz -lboost_program_options -I ./Tools/ --std=c++11

run:
	./bin/parser -i ./src/test/test_data/inFastqTest.fastq -o ./output.csv -p [NNNN][ATCAGTCAACAGATAAGCGA][NNNN] -m 1,2,1 -t 1 -b ./src/test/test_data/barcodeFile.txt
