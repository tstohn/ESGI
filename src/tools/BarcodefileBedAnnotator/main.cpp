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

#include "BarcodeBedAnnotator.hpp"

/** 
 * @brief Simple tool to annotate a TSV file with tab seperated barcodes with a column from a bed file. 
 * E.g.: We have a barcode file from ezgi, the first column MUST contain the read name, this tool then takes an
 * additional bed file and a column index (0-indexed) and adds this column element to the barcode file if the read-name
 * exists in the bed file. (Due to overlapping read-mappings the first occuring gene-name is assigned to the barcode file).
 **/

using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, std::string& barcodeFile, std::string& bedFile, int& featureCol)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input-file,i", value<std::string>(&(barcodeFile))->required(), "Barcodes file, contains all barcodes that were mapped for a fastq-read name.")
            ("bed-file,b", value<std::string>(&(bedFile))->required(), "Bed-file containing read-names (in 4th column) and feature annotations like \
            gene names. The column-index (0-indexed) with the feature must also be given with -f.")
            ("feature-column,f", value<int>(&(featureCol))->required(), "feature column to annotate the barcode file with. E.g., column 14 of the bed file (0-indexed).")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";

            std::cout << "###########################################\n";
            std::cout << "EXAMPLE CALL:\n ";
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
    std::string barcodeFile;
    std::string bedFile;
    int featureCol;
    if(parse_arguments(argv, argc, barcodeFile, bedFile, featureCol))
    {
        BarcodeBedAnnotator barcodeann = BarcodeBedAnnotator(barcodeFile, bedFile, featureCol);
        // run barcode annotation
        barcodeann.annotate();
    }

    return EXIT_SUCCESS;
}
