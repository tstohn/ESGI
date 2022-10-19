#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/version.hpp>

#include "BarcodeProcessingHandler.hpp"
#include "BarcodeMapping.hpp"
#include "MappingAroundLinker.hpp"

/**
 * @brief A Tool to use for quality control of the Combinatorial Indexing run:
 * We map sequentially the linker sequences to a substring which is subread = read[last constant match ... end of read].
 * If the same sequence appear several times in the read, we map ONLY the first occurence (bcs ideally if the same sequence occurs several times, it should
 * also be a repretitive pattern in the ideal sequence.) This means wrongly repetitive sequences are NOT detected.
 * After each mapping of a constant sequence we map sequentially all vairable barcodes which are before this constant barcode. If we do not find 
 * one barcode we continue with the next until we tried all barcodes before the last mapped linker. Those variable barcodes are mapped to a 
 * subsequence = read[end of 2nd last linker match...beginning of last linker match]. After each mapped variable barcode which lies in this sequence, we remove the mapped 
 * bases from the beginning of the substring.
 * 
 */

using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, std::string& inputFile, std::string& output, 
                     int& threads, std::string& barcodeFile,
                     std::string& seqpattern, std::string& mismatches)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("output,o", value<std::string>(&(output))->required(), "output file with all split barcodes")
            ("barcodeList,b", value<std::string>(&(barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC. This can be the same list as it was for FastqParser.")

            //parameters for failed lines
            ("inputFile,i", value<std::string>(&inputFile), "input text file, to analyze barcode pattern.")
            ("sequencePattern,p", value<std::string>(&seqpattern), "pattern for the sequence to match, \
            every substring that should be matched is enclosed with square brackets. N is a wild card match, a substring \
            should not be a combination of wild card chars and constant chars : [AGCTATCACGTAGC][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")
            ("mismatches,m", value<std::string>(&(mismatches)), "list of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2")

            ("thread,t", value<int>(&threads)->default_value(1), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";

            std::cout << "###########################################\n";
            std::cout << "EXAMPLE CALL:\n ./bin/parser -i ./inFile -o ./outFile -p [AGTCAGTC][NNNN] -b ./barcodeFile.txt -m 2,1 -t 5\n";
            std::cout << "Calling the Tool, mapping each line to a pattern, that starts with a constant sequence of <AGTCAGTC> in which up to two mismatches\n\
            are allowed. After that the reads should contain a four base long sequence with a maximum of one mismatch. All sequences that are allowed are in the text file barcodeFile.txt.\n\
            This file should have in its first row, all comma seperated sequences that could map, and they should be all only four bases long, allowed chars are only A,C,G,T.\n";
            std::cout << "###########################################\n";

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

int main(int argc, char** argv)
{
    std::string inputFile;
    std::string output;
    std::string barcodeFile;
    std::string seqpattern;

    int threads;
    std::string mismatches;
    if(parse_arguments(argv, argc, inputFile, output, threads, 
                       barcodeFile, seqpattern, mismatches))
    {

        //analze lines that could not be mapped
        input input;
        input.barcodeFile = barcodeFile;
        input.inFile = inputFile;
        input.outFile = output;
        input.patternLine = seqpattern;
        input.mismatchLine = mismatches;
        input.threads = threads;

        // run demultiplexing
        if(!( endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz") ||  endWith(input.inFile, "txt") ))
        {
            std::cout << "Wrong file format for input file <-i>!\n";
            exit(EXIT_FAILURE);
        }
        
        if(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz"))
        {
            MappingAroundLinker<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromFastqFilePolicy> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "txt"))
        {
            MappingAroundLinker<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy> mapping;
            mapping.run(input);
        }
        else
        {
            fprintf(stderr,"Input file must be of format: <.fastq> | <.fastq.gz> | <.txt>!!!\nFail to open file: %s\n", input.inFile.c_str());
            exit(EXIT_FAILURE);
        }

    }
 
    return EXIT_SUCCESS;
}
