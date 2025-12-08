#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>
#include <string_view>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <regex>
#include <thread>
#include <mutex>
#include <boost/thread/thread.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <cmath>
#include <unordered_map>
#include <filesystem>

#include "seqtk/kseq.h"
#include "dataTypes.hpp"
#include "DemultiplexedLine.hpp"
#include "DemultiplexedStatistics.hpp"
#include "helper.hpp"

KSEQ_INIT(gzFile, gzread)

/** @brief mapping sequentially each barcode leaving no pattern out,
 *if a pattern can not be found the read is discarded
 **/
class MapEachBarcodeSequentiallyPolicy
{
    public:
        bool split_line_into_barcode_patterns(
            const std::pair<fastqLine, fastqLine>& seq, 
            DemultiplexedLine& demultiplexedLine, const input& input,
            BarcodePatternPtr barcodePatterns, 
            int& mmScore, OneLineDemultiplexingStatsPtr stats);
};

/** @brief like the sequential barcode mapping policy, for paired-end reads
 **/
class MapEachBarcodeSequentiallyPolicyPairwise
{
    private:
        bool map_forward(const fastqLine& seq,
                        BarcodePatternPtr barcodePatterns,
                        OneLineDemultiplexingStatsPtr stats,
                        DemultiplexedLine& demultiplexedLine,
                        unsigned int& barcodePosition,
                        int& score_sum,
                        PatternType type);
        bool map_reverse(const fastqLine& seq, 
                        BarcodePatternPtr barcodePatterns,
                        OneLineDemultiplexingStatsPtr stats,
                        DemultiplexedLine& demultiplexedLine,
                        unsigned int& barcodePosition,
                        int& score_sum);
        bool combine_mapping(const BarcodePatternPtr& barcodePatterns,
                             DemultiplexedLine& demultiplexedLineFw, //this list is extended to real list
                             const unsigned int& barcodePositionFw,
                             const DemultiplexedLine& demultiplexedLineRv,
                             const unsigned int& barcodePositionRv,
                             OneLineDemultiplexingStatsPtr stats,
                             OneLineDemultiplexingStatsPtr statsRv,
                             int& score_sum);
    public:
        bool split_line_into_barcode_patterns(
            const std::pair<fastqLine, fastqLine>& seq,  
            DemultiplexedLine& demultiplexedLine,
            const input& input, 
            BarcodePatternPtr barcodePatterns, int& mmScore,
            OneLineDemultiplexingStatsPtr stats);
};

/**
** Mapping Linker (constant) sequences first.
** Linker sequences are mapped to the whole length of the sequence.
**/
class MapAroundConstantBarcodesAsAnchorPolicy
{
    public:
        bool split_line_into_barcode_patterns(
            const std::pair<fastqLine, fastqLine>& seq, 
            DemultiplexedLine& demultiplexedLine,
            const input& input,
            BarcodePatternPtr barcodePatterns, int& mmScore, 
            OneLineDemultiplexingStatsPtr stats);
        void map_pattern_between_linker(const std::string& seq, const int& oldEnd, const int& start, 
                                        BarcodePatternPtr barcodePatterns, std::vector<std::string>& barcodeList,
                                        int& barcodePosition, int& skippedBarcodes);
};

/// parser policy for txt files
class ExtractLinesFromTxtFilesPolicy
{
    public:
    void init_file(const std::string& fwFile, const std::string& rvFile)
    {       
        (void)rvFile; //we have only a forward fastq-read

        // check if file exists
        if (!std::filesystem::exists(fwFile))
        {
            throw std::runtime_error("File does not exist: " + fwFile);
        }
        // check it's a regular file, e.g., not a directory
        if (!std::filesystem::is_regular_file(fwFile))
        {
            throw std::runtime_error("Path is not a regular file: " + fwFile + ". Did you by accident provide a directory?\n");
        }

        //no error handling for txt file right now
        fileStream.open(fwFile, std::ios::in);

        //check if file can be opened & count lines
        if (fileStream.is_open()) 
        {
            fileStream.seekg(0);
            //counting newlines for proxy of read number
            totalReads = std::count(std::istreambuf_iterator<char>(fileStream),
                                    std::istreambuf_iterator<char>(), '\n');
                                        // Check last character (only if file isn't empty)
            //if the last line has no newline add +1
            if (fileStream.tellg() > 0) 
            {
                fileStream.clear(); // Clear EOF flag
                fileStream.seekg(-1, std::ios::end);
                char lastChar;
                fileStream.get(lastChar);
                if (lastChar != '\n') {
                    totalReads += 1;
                }
            }
        } 
        else 
        {
            std::cerr << "Error opening input txt-file!" << std::endl;
            exit(EXIT_FAILURE);
        }

        fileStream.clear();
        fileStream.seekg(0);
    }

    //for txt files we assume every line contains a line of bases
    //quality and read names DO NOT exist
    bool get_next_line(std::pair<fastqLine, fastqLine>& line)
    {   
        bool returnValue = true;
        if(!std::getline(fileStream, line.first.line))
        {
            returnValue = false;
        }
        if(returnValue)
        {
            line.first.line.erase(std::remove(line.first.line.begin(), line.first.line.end(), '\n'), line.first.line.end());
        }

        return(returnValue);
    }

    void close_file()
    {
        fileStream.close();
    }

    unsigned long long get_read_number()
    {
        return totalReads;
    }
    
    std::ifstream fileStream;
    unsigned long long totalReads;
};

///parser policy for fastq(.gz) files
class ExtractLinesFromFastqFilePolicy
{
    public:

    void init_file(const std::string& fwFile, const std::string& rvFile)
    {
        (void)rvFile; //we have only a forward fastq-read

        // check if file exists
        if (!std::filesystem::exists(fwFile))
        {
            throw std::runtime_error("File does not exist: " + fwFile);
        }
        // check it's a regular file, e.g., not a directory
        if (!std::filesystem::is_regular_file(fwFile))
        {
            throw std::runtime_error("Path is not a regular file: " + fwFile + ". Did you by accident provide a directory?\n");
        }

        gzFile fpCount = gzopen(fwFile.c_str(), "r");
        if (fpCount == Z_NULL)
        {
            throw std::runtime_error("Invalid gz file: " + fwFile);
        }
        kseq_t* ksCount = kseq_init(fpCount);

        unsigned long long count = 0;
        while (kseq_read(ksCount) >= 0)
        {
            ++count;
            if (count == ULLONG_MAX)
            {
                std::cerr << "WARNING: too many reads to count accurately; disabling progress bar.\n";
                break;
            }
        }
        kseq_destroy(ksCount);
        gzclose(fpCount);
        totalReads = count;

        // re-open file for real processing
        fp = gzopen(fwFile.c_str(), "r");
        if (fp == Z_NULL)
            throw std::runtime_error("Unable to reopen gz file: " + fwFile);

        ks = kseq_init(fp);
    }

    bool get_next_line(std::pair<fastqLine, fastqLine>& line, bool reverse = false)
    {
        if(kseq_read(ks) < 0)
        {
            if(strlen(ks->seq.s) != strlen(ks->qual.s))
            {
                std::cout << ks->seq.s << "\n";
                std::cout << ks->qual.s << "\n";

                std::cout << "Warning: base quality and read are of different length!\n";
                exit(EXIT_FAILURE);
            }
            return false;
        }
        if(!reverse)
        {
            line.first.line = std::string(ks->seq.s);
            line.first.quality = std::string(ks->qual.s);
            line.first.name = std::string(ks->name.s);
        }
        else
        {
            line.second.line = std::string(ks->seq.s);
            line.second.quality = std::string(ks->qual.s);
            line.second.name = std::string(ks->name.s);
        }
        return true;
    }

    void close_file()
    {
        kseq_destroy(ks);
        gzclose(fp);
    }

    unsigned long long get_read_number()
    {
        return totalReads;
    }

    kseq_t* ks;
    unsigned long long totalReads;
    gzFile fp;

};

class ExtractLinesFromFastqFilePolicyPairedEnd
{

    public:
        void init_file(const std::string& fwFile, const std::string& rvFile)
        {
            fwFileManager.init_file(fwFile, "");
            rvFileManager.init_file(rvFile, "");
        }

        bool get_next_line(std::pair<fastqLine, fastqLine>& line)
        {
            bool fwBool = fwFileManager.get_next_line(line);
            bool rvBool = rvFileManager.get_next_line(line, true);
            return(fwBool&&rvBool);
        }

        void close_file()
        {
            fwFileManager.close_file();
            rvFileManager.close_file();
        }

        unsigned long long get_read_number()
        {
            assert(fwFileManager.get_read_number() == rvFileManager.get_read_number());
            return(fwFileManager.get_read_number());
        }

        ExtractLinesFromFastqFilePolicy fwFileManager;
        ExtractLinesFromFastqFilePolicy rvFileManager;
};

/** @brief generic class for the barcode mapping
 * @param MappingPolicy: the policy used to map one barcode after the other, probably mostly used one should be
 * MapEachBarcodeSequentiallyPolicy
 * @param FilePolicy: the policy used for file reading, txt or fastq(.gz) files
 * @details usage: call the 'run' method, to perform the barcode mapping, calling then get_demultiplexed_reads returns
 * all the mapped barcodes per read as a BarcodeMappingVector, careful, this data is only valid as long as the Mapping is
 **/
template<typename MappingPolicy, typename FilePolicy>
class Mapping : protected MappingPolicy, protected FilePolicy
{
    public:
        ///explicit costructor to initiate the uniqueCharSet for the barcodes in DemultiplexedReads
        Mapping()
        {

            //lock to write progress in terminal
            printProgressLock = std::make_unique<std::mutex>();
            barcodePatternList = std::make_shared<std::vector<BarcodePatternPtr>>();
        }

    private:

        //parses all input arguments from barcode order, to mismatches in which barcodes etc.
        void parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, 
                                std::vector<std::vector<std::string> >& varyingBarcodes);

        //lock for update bar
        std::unique_ptr<std::mutex> printProgressLock;  

    protected:

        //list of all the possible barcode patterns
        MultipleBarcodePatternVectorPtr barcodePatternList;

        //dictionary mapping the barcodePattern-name to its list of succesfully demultiplexed reads
        std::vector<DemultiplexedReads> demultiplexReadsList;

        //process all the input information and check for validity
        //e.g. delete old output files if present, parse barcode from barcodeFile, match them to their
        //number of mismatches etc.
        void initialize_mapping(const input& input);
        
        BarcodePatternPtr create_barcodeVector_from_patternLine(
            const std::vector<std::string>& barcodeList, 
            const std::vector<int>& mismatchList, 
            const std::string& patternName,
            std::unordered_map<std::string, std::vector<std::string>>& fileToBarcodesMap,
            const input& input);
        bool generate_barcode_patterns(const input& input);

//OLD FUNCTION:
        //generate the structure of all reads, which barcode has to be mapped where with how many mismatches
        //basically it is a vector of Barcode objects, this function calls 'parse_barcode_data' and return a vector of
        //pairs that hold <barcode-regex, char determining the kind of barcode> with kind of barcode beeing e.g. a variable, constant, etc.
//  std::vector<std::pair<std::string, char> > generate_barcode_patterns(const input& input);

        //return the structure holder our barcode pattern, that we try to map to every read
        const MultipleBarcodePatternVectorPtr get_barcode_pattern()
        {
            return barcodePatternList;
        }

        //wrapper to call the actual mapping function on one read and updates the status bar
        bool demultiplex_read(const std::pair<fastqLine, fastqLine>& seq, 
                              DemultiplexedLine& demultiplexedLine,
                              BarcodePatternPtr pattern,
                              const input& input, 
                              std::atomic<unsigned long long>& count, const unsigned long long& totalReadCount,
                              int& mmScore,
                              OneLineDemultiplexingStatsPtr stats);

};

//functions to parse pattern file and mismatches - also used in final ESGI tool therefore made accessible here
std::vector<std::pair<std::string, std::vector<std::string>>> parse_pattern_file(const std::string& patternFile);
std::vector<std::vector<int>> parse_mismatch_file(const std::string& mismatchFile);
