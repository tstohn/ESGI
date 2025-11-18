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

        if (vm.count("help")) 
        {
            // 1.) PRINT THE USAGE
            std::cout << "Usage: " << argv[0] << " <input.ini>\n";
            std::cout << desc << "\n";

            // 2.) PRINT AN INI EXAMPLE
            std::filesystem::path exePath = std::filesystem::canonical(argv[0]).parent_path().parent_path(); //the path one above bin
            std::filesystem::path iniPath = exePath / "src/test/test_data/test_esgi/esgi_example.ini";

            std::ifstream in;
            in.open(iniPath);
            if (in) 
            {
                std::cout << "You have to run esgi with a single argument, the input.ini which contains all the information of how to run esgi.\n";
                std::cout << "You find an example below.\n";
                std::cout << 
                "╔═════════════════════╗\n"
                "║   EXAMPLE INI FILE  ║\n"
                "╚═════════════════════╝\n";
            }
            if (!in) 
            {
                std::cerr << "Could not find the example.ini file for esgi.\n";
                return 0;
            }
            std::ostringstream buf;
            buf << in.rdbuf();
            std::cout << buf.str() << "\n";
            
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
    bool dnaPatternPresent = containsDNAPattern(config.patternFile);

    // if we r on windows and the pattern contain RNA throw an error (therefore we need to quickly parse the pattern file)
    if(windows && dnaPatternPresent)
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
    if(!run_demultiplex(config))
    {
        std::cerr << "EXITING ESGI: demultiplexing failed!\n";
    }
    //output file of demultipelxed reads is DIR/PREFIX_(PATTERN_0 | PROTEIN).tsv
    std::vector<std::pair<std::string, std::vector<std::string>>> patterns = parse_pattern_file(config.patternFile);
    std::vector<std::vector<int>> mismatches = parse_mismatch_file(config.mismatchesFile);
    if(patterns.size() > 1 || mismatches.size() > 1){std::cerr << "For now ESGI can handle only a single pattern. If you want to demultiplex several patterns at once\n"\
    "we recommend to use the tools demultiplex and count individually since every pattern might have different positions for different features.";}
    if(patterns.size() != mismatches.size() ){std::cerr << "The number of patterns in pattern file is not equal to the number of mismatches in mismatch file. Even the \
    [-],[*] patterns need an arbitrary mismatch number like 0\n";}
    // if we have no RNA data the extension is tsv, if we do have a RNA pattern we have a tsv and a fastq file
    std::string demultiplexOutput = stripQuotes(patterns.at(0).first);
    if(config.prefix != "")
    {
        demultiplexOutput = config.prefix.value() + "_" + demultiplexOutput;
    }
    //the output of deultiplexing wihtout the file ending: it can be .tsv for patterns or .fastq for DNA
    std::filesystem::path demultiplexOutputPath = std::filesystem::path(config.output) / demultiplexOutput;
    intermediateFiles.demultiplexingOutput = demultiplexOutputPath.string();
    
    // #########################
    // call rna mapping
    // #########################
    if(dnaPatternPresent)
    {
        #ifdef ENABLE_HTSLIB
            std::cout <<
            "╔═════════════════════════════════════════════════════════════════════════════════╗\n"
            "║ 2a.) RUN RNA-MAPPING: run STAR and annotate reads from step 1 with mapped genes ║\n"
            "╚═════════════════════════════════════════════════════════════════════════════════╝\n";
            if(!run_rna_mapping(config, intermediateFiles))
            {
                std::cerr << "EXITING ESGI: running RNA mapping failed!\n";
            }
        #else
            std::cout << "HTSLIB is not installed and needed for the annotation of the STAR output\n";
            std::cout << "HTSLIB does not come with the downloadable binaries of ESGI and has to be installed manually\n";
            std::cout << "Please install htslib or run without a DNA/RNA sequence.\n";
            exit(EXIT_FAILURE);
        #endif
    }

    // #########################
    //call count
    // #########################
    //return a position where we have -,*. Then from every index beyond this we need to substract -1, and no index can be equal to this (e.g., feature, single-cell, UMI indices)
    int specialPatternPos = get_special_pattern_pos(patterns.at(0).second);
    intermediateFiles.specialPatternPos = specialPatternPos;
    int dnaPatternPos = get_DNA_pattern_pos(patterns.at(0).second);
    intermediateFiles.dnaPatternPos = dnaPatternPos;

    intermediateFiles.patternLength = static_cast<int>(patterns.at(0).second.size());
           
    //set the number of MM for the UMI
    if(config.UMI_ID.has_value()){intermediateFiles.umiMismatches = mismatches.at(0).at(std::stoi(config.UMI_ID.value()));}

    // output file for count is basically just the output file of demultiplex, count will then add a COUNTDATA prefix
    std::string demultiplexOutputPatternFile = demultiplexOutput + ".tsv";

    // the output for counting should be the same output name from demultiplexing: count will then add a COUNTDATA_ prefix
    std::filesystem::path countingOutputPath = std::filesystem::path(config.output) / demultiplexOutputPatternFile;
    intermediateFiles.countingOutput = countingOutputPath.string();

    //enable DNA/RNA mapping with STAR and subsequent annotation only if htslib is installed
    if (dnaPatternPresent) {
        std::cout <<
        "╔══════════════════════════════════════════════════════╗\n"
        "║ 2b.) RUN COUNTING: create single-cell feature matrix ║\n"
        "╚══════════════════════════════════════════════════════╝\n";

        // now counting input is the annotate data, this is the origional demultiplexed output with annotated-suffix
        std::string demultiplexOutputAnnotated = demultiplexOutput + "_annotated.tsv";
        std::filesystem::path annotatedOutputPath = std::filesystem::path(config.output) / demultiplexOutputAnnotated;
        intermediateFiles.countingInput = annotatedOutputPath.string();
    }
    else {
        std::cout <<
        "╔═════════════════════════════════════════════════════╗\n"
        "║ 2.) RUN COUNTING: create single-cell feature matrix ║\n"
        "╚═════════════════════════════════════════════════════╝\n";
        //there was no RNA-mapping before, the input to counting is just the output of demultiplexing
        intermediateFiles.countingInput = countingOutputPath.string();
    }
    
    if(!run_count(config, intermediateFiles, dnaPatternPresent))
    {
        std::cerr << "EXITING ESGI: counting failed!\n";
    }

    return EXIT_SUCCESS;
}