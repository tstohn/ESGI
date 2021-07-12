#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>
#include "Barcode.hpp"

/** @param:
 *  sequence pattern: could be a sequence of constant and variable nucleotides: e.g.: [xxxxx][AGCGTACGCGAGT][xxxxx][AAGCGtAGCTTC][xxxxx] 
 *  mismatches per sequence in bracket: e.g.: 1,2,1,2,1
 * 
 * 
 **/

//TODO:
/*
make a basetype of barcode:
 - derived classes of constant codes with ONE string and variables one with X strings
 - iterate through barcodes and match them in increasing order

 - if everything matches within the order of mismatches reporta vector of matches

 INPUT: pattern plus mismatches and barcode sequences
 => vector of barcode pointer which inherits to constant and variable (difference std::stirng anf vector of strings)

 OUTPUT: one vector of (vectors: all matches sequences in a row)

PLAN:
- design classes plus a function to calc dist for one or x sequences and return best match
- parse pattern, mismatches, and barcodeList and store in new Baseclasspointer vector


*/

#include "seqtk/kseq.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

KSEQ_INIT(gzFile, gzread)

//old data types
typedef std::vector<std::string> BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;
using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, input& input)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&(input.inFile))->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&(input.outFile))->required(), "output file with all split barcodes")
            
            ("sequencePattern,p", value<std::string>(&(input.patternLine))->required(), "pattern for the sequence to match, \
            every substring that should be matched is enclosed with square brackets. N is a wild card match, a substring \
            should not be a combination of wild card chars and constant chars : [AGCTATCACGTAGC][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")
            ("barcodeList,b", value<std::string>(&(input.barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC")
            ("mismatches,m", value<std::string>(&(input.mismatchLine))->default_value("1"), "list of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2")
           
            ("threat,t", value<int>(&(input.threads))->default_value(5), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/parser -i ./inFile -o ./outFile -p [AGTCAGTC][NNNN] -m 2,1 -t 5\n";
            return false;
        }

        notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    return true;

}

//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
bool split_line_into_barcode_mappings(const std::string seq, input* input, BarcodeMapping& barcodeMap, fastqStats& stats)
{
    //iterate over BarcodeMappingVector
    //for every barcodeMapping element find a match 

    //add this match to the BarcodeMapping

    return true;
}

void write_file(std::string output, BarcodeMappingVector barcodes)
{
    std::ofstream outputFile;
    outputFile.open (output);
    //write header line
    outputFile << "header\n";
    outputFile.close();
}

void map_pattern_to_fastq_lines(std::vector<std::string> fastqLines, input* input, BarcodeMappingVector& barcodes, fastqStats& stats)
{
    BarcodeMapping barcode;
    for(const std::string& line : fastqLines)
    {
        if(split_line_into_barcode_mappings(line, input, barcode, stats))
        {
            barcodes.push_back(barcode);
        }
    }
}

BarcodeMappingVector iterate_over_fastq(input& input)
{
    //read all fastq lines into str vector
    gzFile fp;
    kseq_t *ks;
    fp = gzopen(input.inFile.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", input.inFile.c_str());
    }
    ks = kseq_init(fp);
    std::vector<std::string> fastqLines;
    while( kseq_read(ks) >= 0 ){
        fastqLines.push_back(std::string(ks->seq.s));
    }

    //split input lines into thread buckets
    std::vector<std::vector<std::string> > fastqLinesVector; //vector holding all the buckets of fastq lines to be analyzed by each thread
                                                            // each vector element is one bucket for one thread
    int element_number = fastqLines.size() / input.threads;
    std::vector<std::string>::iterator begin = fastqLines.begin();
    for(int i = 0; i < input.threads; ++i)
    {
        std::vector<std::string>::iterator end = ( i==input.threads-1 ? fastqLines.end() : begin + element_number);

        std::vector<std::string> tmpFastqLines(begin, end);
        fastqLinesVector.push_back(tmpFastqLines);
        begin = end;
    }

    //for every bucket call a thread
    //tmp variables to store thread results
    std::vector<BarcodeMappingVector> barcodesThreadList(input.threads); // vector of mappedSequences to be filled
    std::vector<fastqStats> statsThreadList(input.threads); // vector of stats to be filled
    std::vector<std::thread> workers;
    //final data variables
    BarcodeMappingVector barcodeVectorFinal;
    fastqStats fastqStatsFinal;
    for (int i = 0; i < input.threads; ++i) {
        workers.push_back(std::thread(map_pattern_to_fastq_lines, fastqLinesVector.at(i), &input, std::ref(barcodesThreadList.at(i)), std::ref(statsThreadList.at(i))));
    }
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    //combine thread data
    for (int i = 0; i < input.threads; ++i) 
    {
        barcodeVectorFinal.insert(barcodeVectorFinal.end(), barcodesThreadList.at(i).begin(), barcodesThreadList.at(i).end());
        fastqStatsFinal.perfectMatches += statsThreadList.at(i).perfectMatches;
        fastqStatsFinal.noMatches += statsThreadList.at(i).noMatches;
        fastqStatsFinal.moderateMatches += statsThreadList.at(i).moderateMatches;
    }
    std::cout << "MATCHED: " << fastqStatsFinal.perfectMatches << " | MODERATE MATCH: " << fastqStatsFinal.moderateMatches
              << " | MISMATCHED: " << fastqStatsFinal.noMatches << "\n";
    kseq_destroy(ks);
    gzclose(fp);
    return barcodeVectorFinal;
}

BarcodePatternVectorPtr generate_barcode_patterns(input input)
{
    BarcodePatternVectorPtr barcodePatternVector;

    std::vector<std::string> patterns; // vector of all string patterns
    std::vector<int> mismatches; // vector of all string patterns

    try{
        // parse the pattern, mismatches, and barcode file (perform quality check as well)
        //PARSE PATTERN
        std::string pattern = input.patternLine;
        std::string delimiter = "]";
        const char delimiter2 = '[';
        size_t pos = 0;
        std::string seq;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            if(seq.at(0) != delimiter2)
            {
                std::cerr << "PARAMETER ERROR: Wrong barcode patter parameter, it looks like a \'[\' is missing\n";
                exit(1);
            }
            seq.erase(0, 1);
            patterns.push_back(seq);
        }
        //PARSE mismatches
        pattern = input.mismatchLine;
        delimiter = ",";
        pos = 0;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            mismatches.push_back(stoi(seq));
        }
        mismatches.push_back(stoi(pattern));
        if(patterns.size() != mismatches.size())
        {
            std::cerr << "PARAMETER ERROR: Number of barcode patterns and mismatches is not equal\n";
            exit(1);
        }
        //PARSE barcode file
        
    }
    catch(std::exception& e)
    {
        std::cerr << "PARAMETER ERROR: Check the input parameters again (barcode pattern list, mismatch list, file with barcode lists)\n";
        std::cerr << e.what() << std::endl;
        exit(1);
    }


    for(int i=0; i < patterns.size(); ++i)
    {
        //std::cout << patterns.at(i) << mismatches.at(i)<< " \n";
    }

    return barcodePatternVector;

}


void split_barcodes(input input, BarcodeMappingVector& barcodes)
{
    BarcodePatternVectorPtr barcodePatterns = generate_barcode_patterns(input);
    //if directory iterate over all fastqs
    barcodes = iterate_over_fastq(input);

    //add plate number and combine all barcodes

}

int main(int argc, char** argv)
{
    input input;
    BarcodeMappingVector barcodes;
    if(parse_arguments(argv, argc, input))
    {
        generate_barcode_patterns(input);
        //split_barcodes(input, barcodes);
        //write_file(input.outFile, barcodes);
    }
 
    return EXIT_SUCCESS;
}