// THIS FILE CONTAINS ALL THE TOOL-CALL NEEDED FOR ESGI:
// like demultiplexing, counting, etc

#include <cstdlib>
#include <string>
#include <iostream>
#include <filesystem>

inline bool call_demultiplex(const ESGIConfig& config)
{

    std::filesystem::path toolPath = std::filesystem::path("bin") / "demultiplex";

    // Build command string with parameters
    std::string cmd = toolPath.string() +
        " --input " + config.forward +
        " --output " + config.output +
        " --threads " + std::to_string(config.threads) +
        " --barcodePatternsFile " + config.patternFile +
        " --mismatchFile " + config.mismatchesFile;

    // arguments with default values
    std::string independentStr = config.independent ? "true" : "false";
    cmd += " --independent " + independentStr;
    std::string hammingStr = config.hamming ? "true" : "false";
    cmd += " --hamming " + hammingStr;
    std::string writeFailedLinesStr = config.writeFailedLines ? "true" : "false";
    cmd += " --writeFailedLines " + writeFailedLinesStr;
    std::string writeStatsStr = config.writeStats ? "true" : "false";
    cmd += " --writeStats " + writeStatsStr;

    // optional arguments
    if(config.reverse.has_value()){cmd += " --reverse " + config.reverse.value();}
    if(config.prefix.has_value()){cmd += " --namePrefix " + config.prefix.value();}
    if(config.fastqReadBucketSize.has_value()){cmd += " --fastqReadBucketSize " + std::to_string(config.fastqReadBucketSize.value());}

    std::cout << "Running demultiplexing as: " << cmd.c_str() << "\n";
    int ret = std::system(cmd.c_str());
    if (ret != 0) 
    {
        std::cout << "DEMULTIPLEX failed with code " << ret << "\n";
        return false;
    }

    return true;
}

inline bool run_star()
{

    return true;
}

inline bool run_annotate()
{
    return true;
}

inline bool run_rna_mapping(const ESGIConfig& config)
{
    //call STAR

    //call annotate

    return true;
}

inline bool call_count(const ESGIConfig& config)
{

    return true;
}