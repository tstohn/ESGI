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
                     int& abIdx, int& treatmentIdx, int& threads)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("inputFails,f", value<std::string>(&inputFails)->required(), "lines of nucleotide sequences that failed the analysis")
            ("inputBarcodes,b", value<std::string>(&inputBarcodeMapping)->required(), "lines of tab seperated barcode mappings")
            ("output,o", value<std::string>(&(output))->required(), "output file with all split barcodes")
            ("antibodyIndex,x", value<int>(&abIdx), "Index used for antibody distinction.")
            ("GroupingIndex,y", value<int>(&treatmentIdx), "Index used to group cells(e.g. by treatment). This is the x-th barcode from the barcodeFile (0 indexed).")
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

void analyse_failed_lines()
{

}

void analyse_same_umis(std::string& inputBarcodeMapping, std::string& output, int& abIdx, int& treatmentIdx, const int& threads)
{
    std::string barcodeFile;
    std::string barcodeIndices;
    //data for protein(ab) and treatment information
    std::string abFile; 
    std::string treatmentFile;
    std::vector<std::string> abBarcodes;
    std::vector<std::string> treatmentBarcodes;
    
    //generate the dictionary of barcode alternatives to idx
    CIBarcode barcodeIdData;
    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData, abBarcodes, abIdx, treatmentBarcodes, treatmentIdx);
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
    int abIdx;
    int treatmentIdx;
    int threads;
    if(parse_arguments(argv, argc, inputFails, inputBarcodeMapping, output, abIdx, treatmentIdx, threads))
    {
        analyse_same_umis(inputBarcodeMapping, output, abIdx, treatmentIdx, threads);
        analyse_failed_lines();
    }
 
    return EXIT_SUCCESS;
}