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

#include "BarcodeBamAnnotator.hpp"

/** 
 * @brief Simple tool to annotate a TSV file with tab seperated barcodes with a column from a BAM file. 
 * E.g.: We have a barcode file from ezgi, the first column MUST contain the read name, this tool then takes an
 * additional BAM file (e.g., from STAR annotations) and a tag that should be used as feature (e.g., GX for gene ids, or GN for gene names) and adds this column element to the barcode file if the read-name
 * exists in the BAM file.
 * Output are only reads that: 1.) Could be mapped the the barcode pattern, 2.) mapped to a gene in the reference genome
 **/

using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, std::string& barcodeFile, std::string& bamFile, std::string& featureTag)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input-file,i", value<std::string>(&barcodeFile)->required(), "Barcodes file, contains all barcodes that were mapped for a fastq-read name.")
            ("bam-file,b", value<std::string>(&bamFile)->required(), "Bam-file with read-mapping. E.g., from running STAR, it must contain the\
            read-namke in the first column (to map the read back to the barcode-reads. And it must contain the feature-tag, given with the -f argument, E.g., GX for gene ids.)")
            ("feature-tag,f", value<std::string>(&featureTag)->default_value("GX"), "feature-tag to annotate the barcode file with. Default are gene ids (GX), for gene names set it to \'GN\'.")

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
    std::string bamFile;
    std::string featureTag;
    if(parse_arguments(argv, argc, barcodeFile, bamFile, featureTag))
    {
        BarcodeBamAnnotator barcodeann = BarcodeBamAnnotator(barcodeFile, bamFile.c_str(), featureTag.c_str());
        // run barcode annotation
        barcodeann.annotate();
    }

    return EXIT_SUCCESS;
}
