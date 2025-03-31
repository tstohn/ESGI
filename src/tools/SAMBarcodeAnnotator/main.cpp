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

/** 
 * @brief A tool to map a demultiplexed file with read_names and barcodes to an annotated BAM file
 **/

using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, std::string& barcodeFile, std::string& bedFile, std::string& featureName)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input-file,i", value<std::string>(&(barcodeFile))->required(), "Barcodes file, contains all barcodes that were mapped for a fastq-read name.")
            ("bed-file,b", value<std::string>(&(bedFile))->required(), "Bed-file containing read-names (in 4th column) and feature annotations like \
            gene names. The column with the feature name must have the feature in quotes: e.g. <featureName \"SERT\">.")

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
    std::string featureName;
    if(parse_arguments(argv, argc, barcodeFile, bedFile, featureName))
    {
        barcodeann = BarcodeBedAnnotator(barcodeFile, bedFile, featureName);
        // run barcode annotation
        barcodeann.annotate();
    }

    return EXIT_SUCCESS;
}
