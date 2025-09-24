#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>

#include <esgiUtils.hpp>
#include <parsingUtils.hpp>

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
    
    //PARSE PARAMETERS
    std::string input;
    if (!parse_arguments(argc, argv, input)) 
    {
        return EXIT_FAILURE;
    }
    if (!validate_arguments(input)) 
    {
        return EXIT_FAILURE;
    }

    //parse the ini file and store all necessary parameter in structure

    // if we r on windows and the pattern contain RNA throw an error

    std::cout << "Input file (" << input << ") is valid - starting esgi.\n";

    // call demultiplex

    //call STAR

    //call annotate

    //call count


    return EXIT_SUCCESS;
}