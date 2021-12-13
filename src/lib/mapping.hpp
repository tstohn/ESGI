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

#include "Barcode.hpp"
#include "seqtk/kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef std::vector< std::shared_ptr<std::string> > BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;

//mapping sequentially each barcode leaving nmo pattern out,
//if a pattern can not be found the read is idscarded
class MapEachBarcodeSequentiallyPolicy
{
    private:
        bool check_if_seq_too_short(const int& offset, const std::string& seq);
    public:
        bool split_line_into_barcode_patterns(const std::string& seq, const input& input, BarcodeMapping& barcodeMap, BarcodeMapping& realBarcodeMap,
                                      fastqStats& stats,
                                      BarcodePatternVectorPtr barcodePatterns);
};

//mapping only constant barocdes as anchor first
class MapAroundConstantBarcodesAsAnchorPolicy
{
    public:
    bool split_line_into_barcode_patterns(const std::string& seq, const input& input, BarcodeMapping& barcodeMap, BarcodeMapping& realBarcodeMap,
                                      fastqStats& stats,
                                      BarcodePatternVectorPtr barcodePatterns);
};

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

//class for barcode mappings, the mapping algorithm is chosen by the 
//mapping policy
template<typename MappingPolicy, typename FilePolicy>
class Mapping : private MappingPolicy, private FilePolicy
{
    private:

        void parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, 
                                std::vector<std::vector<std::string> >& varyingBarcodes);


        BarcodeMappingVector sequenceBarcodes; // uncorrectedSequenceBarcode for later
        BarcodeMappingVector realBarcodes; // vector or all the patterns that we map
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
        void demultiplex_read(const std::string& seq, const input input);
        //run the actual mapping by using MappingPolicy
        void run_mapping(const input& input);


    public:
        void run(const input& input);
};
