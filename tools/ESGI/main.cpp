#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>

#include <ESGIConfig.hpp>
#include <ESGIHelper.hpp>

/**
 * @brief ESGI combines the individual tools in this library, namely demultiplex/ (run STAR with annotate)/ and count
 * ESGI demultiplexes any pattern in fastq files and then creates a single-cell feature matrix.
 * Input is a file.ini, this file MUST have the .ini ending and contains all the important information to demultiplex the data
 * an example file can be found in ./src/test/test_data/test_esgi/esgi_example.ini and additional files in the same directory
**/

namespace fs = std::filesystem;
namespace po = boost::program_options;

// Parse command-line arguments
bool parse_arguments(int argc, char** argv, std::string& input) {
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("input", po::value<std::string>(&input), "input .ini file (positional)");

        po::positional_options_description p;
        p.add("input", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .positional(p)
                      .run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: " << argv[0] << " <input.ini>\n";
            std::cout << desc << "\n";
            return false;
        }

        po::notify(vm);

    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << "\n";
        return false;
    }

    return true;
}

// Validate arguments (file existence and extension)
bool validate_arguments(const std::string& input) 
{
    if (input.empty()) {
        std::cerr << "Error: input .ini file required: this file must end with <.ini> and contains all the information about the experiment that should be demultiplexed. \n";
        std::cerr << "RUN THE TOOL WITH: \"esgi input.ini\" \n";
        std::cerr << "FOR HELP RUN: \"esgi --help\" \n";

        return false;
    }

    if (input.size() < 4 || input.substr(input.size() - 4) != ".ini") 
    {
        std::cout << "Warning: input file does not have .ini extension - as our examples have.\n";
        std::cout << "It is not a requirement as long as you make sure the content of the file follows our example.\n";
    }

    if (!fs::exists(input)) 
    {
        std::cerr << "Error: the parameter (.ini) file does not exist: " << input << "\n";
        return false;
    }

    return true;
}


int main(int argc, char** argv) 
{
    bool windows = false;
    #if defined(_WIN32)
    windows = true;
    #endif

    //PARSE PARAMETERS
    std::string input;
    if (!parse_arguments(argc, argv, input)) 
    {
        return EXIT_FAILURE;
    }
    //check validity of input file
    if (!validate_arguments(input)) 
    {
        return EXIT_FAILURE;
    }

    //parse the ini file and store all necessary parameter in structure
    ESGIConfig config = loadESGIConfigFromFile(argv[1]);
    config.write();
    bool rnaPatternPresent = containsRNAPattern(config.patternFile);

    // if we r on windows and the pattern contain RNA throw an error (therefore we need to quickly parse the pattern file)
    if(windows && rnaPatternPresent)
    {
        std::cerr << "We do not support demultiplexing with a RNA pattern on Windows for now, \n" \
        "We struggled to use the htslib for annotating demultiplexed reads with the STAR output under windows and disabled it for now. \n" \
        "If there is mayor interest in the feature please reach out to us so we can have another look at it.";
    }
    std::cout << "Input file (" << input << ") is valid - starting ESGI.\n";

    // tools like demultiplex, STAR, etc. create intermediate files that are then input for subsequent tools
    //this structure keeps track of those
    ESGIIntermediateFiles intermediateFiles;

    // #########################
    // call demultiplex
    // #########################
    std::cout <<
    "╔══════════════════════════════════════════════════════╗\n"
    "║ 1.) RUN DEMULTIPLEXING: split FASTQ reads into chunks║\n"
    "╚══════════════════════════════════════════════════════╝\n";
    if(!call_demultiplex(config))
    {
        std::cerr << "EXITING ESGI: demultiplexing failed!\n";
    }
    //output file of demultipelxed reads is DIR/PREFIX_(PATTERN_0 | PROTEIN).tsv

    // #########################
    // call rna mapping
    // #########################
    if(rnaPatternPresent)
    {
        std::cout <<
        "╔═════════════════════════════════════════════════════════════════════════════════╗\n"
        "║ 2a.) RUN RNA-MAPPING: run STAR and annotate reads from step 1 with mapped genes ║\n"
        "╚═════════════════════════════════════════════════════════════════════════════════╝\n";
        if(!run_rna_mapping(config))
        {
            std::cerr << "EXITING ESGI: running RNA mapping failed!\n";
        }
    }

    // #########################
    //call count
    // #########################
    if (rnaPatternPresent) {
        std::cout <<
        "╔══════════════════════════════════════════════════════╗\n"
        "║ 2b.) RUN COUNTING: create single-cell feature matrix ║\n"
        "╚══════════════════════════════════════════════════════╝\n";
    }
    else {
        std::cout <<
        "╔═════════════════════════════════════════════════════╗\n"
        "║ 2.) RUN COUNTING: create single-cell feature matrix ║\n"
        "╚═════════════════════════════════════════════════════╝\n";
    }
    if(!call_count(config))
    {
        std::cerr << "EXITING ESGI: counting failed!\n";
    }

    return EXIT_SUCCESS;
}