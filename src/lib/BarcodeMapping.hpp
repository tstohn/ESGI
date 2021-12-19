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
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <cmath>

#include "Barcode.hpp"
#include "seqtk/kseq.h"
#include "dataTypes.hpp"

KSEQ_INIT(gzFile, gzread)

typedef std::vector< std::shared_ptr<std::string> > SequenceMapping;
typedef std::vector<const char*> BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;

/** @brief representation of all the mapped barcodes:
 * basically a vector of all reads, where each read itself is a vector of all mapped barcodes
 * This structures stores each barcode only once, handled by the UniqueCharSet, by that
 * most highly redundant datasets can be stored in only a fraction of its origional memory
**/
class DemultiplexedReads
{
    public:

        DemultiplexedReads()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
            lock = std::make_unique<std::mutex>();
        }

        void addVector(std::vector<std::string> barcodeVector)
        {
            std::lock_guard<std::mutex> guard(*lock);
            BarcodeMapping uniqueBarcodeVector;
            for(std::string barcode : barcodeVector)
            {
                uniqueBarcodeVector.emplace_back(uniqueChars->getUniqueChar(barcode.c_str()));
            }
            mappedBarcodes.push_back(uniqueBarcodeVector);
        }

        const size_t size()
        {
            return mappedBarcodes.size();
        }

        const BarcodeMapping at(const int& i)
        {
            return(mappedBarcodes.at(i));
        }

        const BarcodeMappingVector get_all_reads()
        {
            return(mappedBarcodes);
        }

    private:
        BarcodeMappingVector mappedBarcodes;
        //all the string inside this class are stored only once, 
        //set of all the unique barcodes we use, and we only pass pointers to those
        std::shared_ptr<UniqueCharSet> uniqueChars;
        std::unique_ptr<std::mutex> lock;  

};

/** @brief mapping sequentially each barcode leaving no pattern out,
 *if a pattern can not be found the read is discarded
 **/
class MapEachBarcodeSequentiallyPolicy
{
    private:
        bool check_if_seq_too_short(const int& offset, const std::string& seq);
    public:
        bool split_line_into_barcode_patterns(const std::string& seq, const input& input, DemultiplexedReads& barcodeMap,
                                      BarcodePatternVectorPtr barcodePatterns, fastqStats& stats);
};

/**
** Mapping Linker (constant) sequences first.
** Linker sequences are mapped to the whole length of the sequence.
**/
class MapAroundConstantBarcodesAsAnchorPolicy
{
    public:
    bool split_line_into_barcode_patterns(const std::string& seq, const input& input, DemultiplexedReads& barcodeMap,
                                      BarcodePatternVectorPtr barcodePatterns, fastqStats& stats);
};

/// parser policy for txt files
class ExtractLinesFromTxtFilesPolicy
{
    public:
    void init_file(const std::string& inFile)
    {       
        //no error handling for txt file right now
        fileStream.open(inFile);

        totalReads = std::count(std::istreambuf_iterator<char>(fileStream), std::istreambuf_iterator<char>(), '\n');
        fileStream.clear();
        fileStream.seekg(0);
    }

    bool get_next_line(std::string& line)
    {   bool returnValue = true;
        if(!std::getline(fileStream, line))
        {
            returnValue = false;
        }
        if(returnValue)
        {
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
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
    void init_file(const std::string& inFile)
    {
        fp = gzopen(inFile.c_str(),"r");
        if(NULL == fp)
        {
            fprintf(stderr,"Fail to open file: %s\n", inFile.c_str());
            exit(EXIT_FAILURE);
        }
        totalReads = 0;
        unsigned char buffer[1000];
        while(!gzeof(fp))
        {
            gzread(fp, buffer, 999);
            for (const char &c : buffer) 
            {
                if ( c == '\n' )
                {
                    if(totalReads == ULLONG_MAX)
                    {
                        std::cout << "WARNING: Analysing more than " << std::to_string(ULLONG_MAX) << " reads. There will be no status update\n";
                        gzrewind(fp);
                        ks = kseq_init(fp);
                        return;
                    }
                    ++totalReads;
                }
            }
        }
        totalReads = totalReads/4;
        gzrewind(fp);
        ks = kseq_init(fp);
    }

    bool get_next_line(std::string& line)
    {
        if(kseq_read(ks) < 0)
        {
            return false;
        }

        line = std::string(ks->seq.s);
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
            barcodeMap = DemultiplexedReads();
            printProgressLock = std::make_unique<std::mutex>();
            stats.statsLock = std::make_unique<std::mutex>();
        }

        /** @brief returns our mapped reads, its a vector (for all raeds) of a vector (for all barcodes sequentially)
         * of const char*, only valid as long as Mapping object - handle with care!!!!
         **/
        const BarcodeMappingVector get_demultiplexed_reads()
        {
            return(barcodeMap.get_all_reads());
        }
        ///number of perfect matches
        const unsigned long long get_perfect_matches()
        {
            return stats.perfectMatches;
        }
        ///number of matches with mismatches
        const unsigned long long get_moderat_matches()
        {
            return stats.moderateMatches;
        }
        ///number of failed matches
        const unsigned long long get_failed_matches()
        {
            return stats.noMatches;
        }
        ///the dictionary of mismatches per barcode
        const std::map<std::string, std::vector<int> > get_mismatch_dict()
        {
            return stats.mapping_dict;
        }

        ///run mapping over all reads of the input file
        void run(const input& input);

    private:

        //parses all input arguments from barcode order, to mismatches in which barcodes etc.
        void parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, 
                                std::vector<std::vector<std::string> >& varyingBarcodes);

        /** @brief the structure that is filled during mapping, kind of our 'end result', that keeps that of all mapped barcodes
        *basically a vector of a vector of mapped barcodes, we can get this data after 'run' by calling 'get_demultiplexed_reads'
        **/
        DemultiplexedReads barcodeMap;

        //representation of the barcode pattern we want to map to all reads
        //basically a vector of Barcode objects (stores all possible barcodes, mismatches that are allowed, etc.)
        BarcodePatternVectorPtr barcodePatterns;
        //statistics of the mapping
        fastqStats stats;
        //lock for update bar
        std::unique_ptr<std::mutex> printProgressLock;  

    protected:

        //initialize the stats dictionary of mismatches per barcode
        void initializeStats();
        //process all the input information and check for validity
        //e.g. delete old output files if present, parse barcode from barcodeFile, match them to their
        //number of mismatches etc.
        void initialize_mapping(const input& input);

        //return the structure holder our barcode pattern, that we try to map to every read
        const BarcodePatternVectorPtr get_barcode_pattern_vector()
        {
            return barcodePatterns;
        }

        //generate the structure of all reads, which barcode has to be mapped where with how many mismatches
        //basically it is a vector of Barcode objects, this function calls 'parse_barcode_data' and return a vector of
        //pairs that hold <barcode-regex, char determining the kind of barcode> with kind of barcode beeing e.g. a variable, constant, etc.
        std::vector<std::pair<std::string, char> > generate_barcode_patterns(const input& input);
        //wrapper to call the actual mapping function on one read and updates the status bar
        bool demultiplex_read(const std::string& seq, const input& input, std::atomic<unsigned long long>& count, const unsigned long long& totalReadCount);
        //run the actual mapping
        void run_mapping(const input& input);
};
