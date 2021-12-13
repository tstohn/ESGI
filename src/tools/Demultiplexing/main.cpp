#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include "Barcode.hpp"
#include "mapping.hpp"

/** 
 * A tool to map fastq lines (stitched to one read, use e.g. fastq-join) to a certain barcode pattern,
 * this pattern must be given as an input parameter plus additional files for patterns with multiple
 * possible barcodes
 * 
 * HOW TO ALGORITHM:
    iterate over patterns and match it to substring plus/minus mismatches on both sides
    allow mismatches at beginning and end (plus/minus those mismatches, bcs imagine a barcode match with a mismatch in the end,
    we then donnt know if its  really a msimatch, or a deletion of the barcode and already part of the next barcode...)
    each matched barcodes is described by the first and last match of the sequence (therefore can be shorter, than real sequence, but not longer)
    UMI or WildcardBarcodes are matched according to the two last matches in the neighboring sequences
    aligning by semi global alignment, bcs it could be that our pattern matches beyond the sequence, that is checked afterwards
 * @param <input> input fastq file (gzipped), considers only full length fastq reads, FW/RV must be stitched together in advance
 * @param <output> output extension, that will be added to the output files, see return
 * @param <sequencePattern> a string with all the barcode patterns, each pattern is enclosed by suqare brackets, valid chars are AGTC ofr bases, N for a sequences
 *                          that hold a fixed number of different combinations, must be declared in the barcodeList file, and X for wildcard sequences, that can
 *                          be anything (e.g. for UMIs): e.g.: [XXXXXXXXXX][AGCGTACGCGAGT][NNNNNNNN][AAGCGtAGCTTC][NNNNNNNN] 
 * @param <barcodeList> a file of discrete barcodes, that can be mapped to the [N...] sequences, each row is one pattern, they r in the order as in the sequence pattern and
 *                      barcodes must be comma seperated
 * @param <mismatches> comma seperated list of mismatches, one entry for each sequence pattern declared above: e.g.: 1,2,1,2,1. This also requires a parameter
                       for the 'X'-sequence, although this parameter is at the moment unused, as we have no sequence to alig this one to
 *  
 * @return  writes a tab seperated file with barcodes, a file with statistical information, and if flag <-r> set, also a tab seperated file
 *          of corresponding reads with their real sequence, including mismatches.
 *          After running a line of perfect, moderate and mismatches is printed. Each count refers to a whole line (all barcodes matched perfectly,
 *          at least one matched only with mismatches, at least one did not match at all. The number of duplicate barcodes is per barcode, and thus might
 *          be greater than the total number of mismatches.)
 * 
 *          barcode file: Output file starts with 'BarcodeMapping_' and ends with the ending extension declared with the -o paramter
 *          stats file: a list where for each barcode in the fastq file, the number of found mismatches is written, the columns start with zero mismatches
 *          to the maximum number of mismtaches declared in mismatches parameter plus one, this last columns sums up all cases of more than max-num mismatches.
 *          real barcode file: like barcode file with uncorrected sequences, starts with RealBarcodeSequence
 **/

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

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
            ("realSequences,r", value<bool>(&(input.storeRealSequences))->default_value(false), "set flag if next to the mapped barcodes also a table of the actual \
            reads - including mismatches - should be printed\n")
            ("fastqReadBucketSize,s", value<int>(&(input.fastqReadBucketSize))->default_value(10000000), "number of lines of the fastQ file that should be read into RAM \
            and be processed, before the next batch of fastq reads is read into RAM and processed.")
            ("failureAssessment,f", value<bool>(&(input.analyseUnmappedPatterns))->default_value(false), "set flag for a quality assessment: In this case patterns that can not be mapped \
            will simply be skipped, to see if reads might follow a wrong pattern\n")

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
    input input;
    if(parse_arguments(argv, argc, input))
    {
        if(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz"))
        {
            Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "txt"))
        {
            Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy> mapping;
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
