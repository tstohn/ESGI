#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <vector>
#include <thread>
#include <pthread.h>
#include <unordered_map>
#include <sstream>
#include <climits>
#include <mutex>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <cmath>

#include "UmiData.hpp"

struct CIBarcode
{
    std::vector<std::unordered_map<std::string, int> > barcodeIdDict; // map of barcode to id, maps are in order of their occurence in the fastqRead
    // ids of the CI barcode within the barcode file (includes only barcodes for variable sequence regions)
    std::vector<int> ciBarcodeIndices; // more or less only of temporary usage, during generation of barcode map vector
    int tmpTreatmentIdx;
};

//statistics of the UMI occurences
struct StatsUmi
{
    //dict: key = number of mismatches; 
    //value = how often this particular number of mismatches between any two UMIs of same AbSc-idx was observed
    std::unordered_map<int, unsigned long long> umiMismatchDict;
};

struct umiQuality
{
    unsigned long long sameUmiDiffAbSc = 0;
    unsigned long long sameUmiSameAbSc = 0;
};

struct umiQualityExtended
{
    unsigned long long numOfUmiOccurences = 0;
    unsigned long long numOfABDifferences = 0;

    unsigned long long numOfBC4Differences = 0;
    unsigned long long numOfBC1Differences = 0;
    unsigned long long numOfBC2Differences = 0;
    unsigned long long numOfBC3Differences = 0;
};

void generateBarcodeDicts(std::string barcodeFile, std::string barcodeIndices, CIBarcode& barcodeIdData, 
                          std::vector<std::string>& proteinDict, const int& protIdx, std::vector<std::string>* treatmentDict = nullptr, const int& treatmentIdx = 0);

class UmiDataParser
{


    public:

        UmiDataParser(CIBarcode barcodeIdData) : barcodeDict(barcodeIdData){}

        void parseFile(const std::string fileName, const int& thread);

        //correct the umis, so that UmiData holds the same UMI for UMIs within a certain mismatch range
        //concurrently fills AbData, while iterating over AbSc lines, count the Ab occurence if this UMI has not been seen before
        void processBarcodeMapping(const int& umiMismatches, const int& thread);

        void writeStats(std::string output);
        void writeUmiCorrectedData(const std::string& output);
        inline void addTreatmentData(std::unordered_map<std::string, std::shared_ptr<std::string> > map)
        {
            rawData.setTreatmentDict(map);
        }
        inline void addProteinData(std::unordered_map<std::string, std::shared_ptr<std::string> > map)
        {
            rawData.setProteinDict(map);
        }

        void extended_umi_quality_check(const int& thread, const std::string& output);

    private:

        //parse the file, store each line in our data structure
        void addFastqReadToUmiData(const std::string& line, const int& elements);
        void parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads);
        void correctUmis(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<abLine>& abDataTmp, 
                         const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected);
        void correctUmisWithStats(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<abLine>& abDataTmp, 
                         const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected);

        void umiQualityCheck(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, umiQuality& qualTmp, int& currentUmisChecked);
        void umiQualityCheckExtended(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, umiQuality& qualTmp, int& currentUmisChecked, const std::string& output);

        //get positions of CIBarcodes
        void getCiBarcodeInWholeSequence(const std::string& line, int& barcodeElements);
        //map all the barcodes of CI to a unique 'number' string as SingleCellIdx
        std::string generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes);

        //data structure storing lines with: UMI, AB_id, SingleCell_id
        //this is the raw data not UMI corrected
        UmiData rawData;

        //new UMI data beeing filled with correcte umi codes
        std::vector<dataLinePtr> umiData;
        //data structure storing lines with: AB_id, SingleCell_id, Ab_count
        std::vector<abLine> abData;

        //Dictionary used to generate the dataLines, maps for each barcode in the sequence all
        //possibilities to an idx
        CIBarcode barcodeDict;
        //statistics of the whole process
        StatsUmi stats;    
        umiQuality qual;  

        std::mutex lock;  
        
        // variables to read the data from each fastq line
        std::vector<int> fastqReadBarcodeIdx; // ids of the CI barcode within the whole string of all barcodes
        int abIdx = INT_MAX;
        int umiIdx = INT_MAX;
        int treatmentIdx = INT_MAX;
        int umiLength = 0;
};
