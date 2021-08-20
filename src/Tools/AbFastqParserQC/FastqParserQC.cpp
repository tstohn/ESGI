#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "UmiDataParser.hpp"

using namespace boost::program_options;

/*
    1.) UMI check
         if same UMI and NOT same AB-SC: check if same AB, check if same BC 4 1 2 3

    2.) Barcode Pattern Check



*/

bool parse_arguments(char** argv, int argc, std::string& inputFails, std::string& inputBarcodeMapping, std::string& output, 
                     int& abIdx, int& treatmentIdx, int& threads, std::string& barcodeFile, std::string& barcodeIndices,
                     std::string& seqpattern, int& mismatches)
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
            ("inputFails,f", value<std::string>(&inputFails), "lines of nucleotide sequences that failed the analysis")
            ("sequencePattern,p", value<std::string>(&seqpattern), "pattern for the sequence to match, \
            every substring that should be matched is enclosed with square brackets. N is a wild card match, a substring \
            should not be a combination of wild card chars and constant chars : [AGCTATCACGTAGC][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")
            ("mismatches,m", value<int>(&(mismatches))->default_value(1), "list of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2")

            //parameters for UMI CHECK
            ("inputBarcodes,i", value<std::string>(&inputBarcodeMapping)->required(), "lines of tab seperated barcode mappings")
            ("antibodyIndex,x", value<int>(&abIdx), "Index used for antibody distinction.")
            ("CombinatorialIndexingBarcodeIndices,c", value<std::string>(&(barcodeIndices))->required(), "comma seperated list of indexes, that are used during \
            combinatorial indexing and should distinguish a unique cell. Be aware that this is the index of the line inside the barcodeList file (see above). \
            This file ONLY includes lines for the varying sequences (except UMI). Therefore the index is not the same as the position in the whole sequence \
            if constant or UMI-seq are present. Index starts with zero.")

            ("thread,t", value<int>(&threads)->default_value(5), "number of threads")

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

void analyse_failed_lines(const std::string& inputFails)
{

}

void analyse_same_umis(std::string& inputBarcodeMapping, std::string& output, int& abIdx, const int& threads,
                    const std::string& barcodeFile, const std::string& barcodeIndices)
{
    //data for protein(ab) and treatment information
    std::string abFile; 
    std::string treatmentFile;
    std::vector<std::string> abBarcodes;

    //generate the dictionary of barcode alternatives to idx
    CIBarcode barcodeIdData;

    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData, abBarcodes, abIdx);

    UmiDataParser dataParser(barcodeIdData);
    //parse the information

    dataParser.parseFile(inputBarcodeMapping, threads);

    //analyse all UMIS
    dataParser.extended_umi_quality_check(threads, output);
}

int main(int argc, char** argv)
{
    std::string inputFails;
    std::string inputBarcodeMapping;
    std::string output;

    std::string barcodeFile;
    std::string barcodeIndices;
    std::string seqpattern;

    int abIdx;
    int treatmentIdx;
    int threads;
    int mismatches;
    if(parse_arguments(argv, argc, inputFails, inputBarcodeMapping, output, abIdx, treatmentIdx, threads, barcodeFile, barcodeIndices, seqpattern, mismatches))
    {
        analyse_same_umis(inputBarcodeMapping, output, abIdx, threads, barcodeFile, barcodeIndices);
        analyse_failed_lines(inputFails);
    }
 
    return EXIT_SUCCESS;
}