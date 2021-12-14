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

#include "Barcode.hpp"
#include "seqtk/kseq.h"
#include "dataTypes.hpp"

KSEQ_INIT(gzFile, gzread)

typedef std::vector< std::shared_ptr<std::string> > SequenceMapping;

typedef std::vector<const char*> BarcodeMapping;

typedef std::vector<BarcodeMapping> BarcodeMappingVector;


//representation of all the mapped barcodes:
//basically a vector of all lines, that contain all barcodes for a read
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

//mapping sequentially each barcode leaving nmo pattern out,
//if a pattern can not be found the read is idscarded
class MapEachBarcodeSequentiallyPolicy
{
    private:
        bool check_if_seq_too_short(const int& offset, const std::string& seq);
    public:
        bool split_line_into_barcode_patterns(const std::string& seq, const input& input, DemultiplexedReads& barcodeMap,
                                      BarcodePatternVectorPtr barcodePatterns);
};

//mapping only constant barocdes as anchor first
/*
class MapAroundConstantBarcodesAsAnchorPolicy
{
    public:
    bool split_line_into_barcode_patterns(const std::string& seq, const input& input, DemultiplexedReads& barcodeMap,
                                      BarcodePatternVectorPtr barcodePatterns);
};*/

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
    
    std::ifstream fileStream;
    int totalReads;
};

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

    kseq_t* ks;
    int totalReads;
    gzFile fp;

};

//keep track of all the mapping results to write for Demultipelx tool
//stats are constructed after splitting fucntion (fill mismatch object, fill morderate,perfect match object)
struct MappingResult
{
    int score;
    SequenceMapping sequences;
};

//class for barcode mappings, the mapping algorithm is chosen by the 
//mapping policy
template<typename MappingPolicy, typename FilePolicy>
class Mapping : private MappingPolicy, private FilePolicy
{
    public:
    //explicit costructor to initiate the uniqueCharSet for the barcodes
        Mapping()
        {
            barcodeMap = DemultiplexedReads();
        }

        const BarcodeMappingVector get_demultiplexed_reads()
        {
            return(barcodeMap.get_all_reads());
        }

    private:

        void parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, 
                                std::vector<std::vector<std::string> >& varyingBarcodes);


//these two SHALL BECOME OBSOLETE AND BE DELETED
        //BarcodeMappingVector sequenceBarcodes; // uncorrectedSequenceBarcode for later
        BarcodeMappingVector realBarcodes; // vector or all the patterns that we map


        DemultiplexedReads barcodeMap;
        BarcodePatternVectorPtr barcodePatterns; // representation of the pattern structure we use fopr mapping

        std::shared_ptr<fastqStats> fastqStatsPtr;

    protected:

        void map_pattern_to_fastq_lines(std::vector<std::string>& fastqLines, const input& input, BarcodeMappingVector& barcodes, 
                                BarcodeMappingVector& realBarcodes, fastqStats& stats, BarcodePatternVectorPtr barcodePatterns,
                                std::vector<std::string>& failedLines);
        void ditribute_jobs_to_threads(const input& input, std::vector<std::string>& fastqLines);


        //INITIALIZATION FUCNTIONS
        const BarcodePatternVectorPtr get_barcode_pattern_vector()
        {
            return barcodePatterns;
        }
        std::vector<std::pair<std::string, char> > generate_barcode_patterns(const input& input);
        //process all the input information and check for validity
        //e.g. delete old output files if present, parse barcode from barcodeFile, match them to their
        //number of mismatches etc.
        void initialize_mapping(const input& input);

        //BARCODE MAPPING FUNCTIONS
        void demultiplex_read(const std::string& seq, const input& input, std::atomic<int>& count, const int& totalReadCount);
        //run the actual mapping by using MappingPolicy
        void run_mapping(const input& input);

    public:
        void run(const input& input);
};
