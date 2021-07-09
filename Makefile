#DEPENDANCIES: zlib
install:
	#download and compile kseq
	mkdir Tools
	cd Tools
	git clone https://github.com/lh3/seqtk/blob/master/kseq.h --branch v1.3
	cd seqtk
	make
	cd ..


parser:
	g++ src/Tools/FastqParser.cpp -o bin/parser -lz -lboost_program_options -I ./Tools/ --std=c++11

run:
	./bin/parser -i ../../TRIALS/Saturation/Plate3SubsamplingKlaas/Subsampled/10.fastq -o ./output.csv -a ATCAGTCAACAGATAAGCGA -w w --antibody ab -t 5