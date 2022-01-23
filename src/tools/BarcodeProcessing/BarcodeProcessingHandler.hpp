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
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <cmath>

#include "DemultiplexedData.hpp"


/**
 * @brief Structure storing a vector with a mapping of the barcode-sequence to a unique ID
 *        for this sequence (this is done seperately for each barcoding round)
 * 
 *        It also stores also the indices for Ci-barcodes, AB, treatment.
 *        These indices indicate the position within the file of all barcodes. So only the position
 *        among all barcodes with the [NNN...] pattern.
 */
struct NBarcodeInformation
{
    // map of barcode to id, maps are in order of their occurence in the fastqRead
    // ids of the CI barcode within the barcode file (includes only barcodes for variable sequence regions)
    std::vector<std::unordered_map<std::string, int> > barcodeIdDict; //used to generate a SC-index
    
    //different indices, they r the index of the NNN-barcodes only (and therefore differe from the index of all barcodes incl. constant ones)
    std::vector<int> NBarcodeIndices; //CI barcode indices
    int NTreatmentIdx; // treatment index
    int NAbIdx; //AB index
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

    std::vector<unsigned long long> AbBarcodeDistribution;
    std::vector<unsigned long long> BC4Distribution;
    std::vector<unsigned long long> BC1Distribution;
    std::vector<unsigned long long> BC2Distribution;
    std::vector<unsigned long long> BC3Distribution;
};

/**
 * @brief Map barcode sequences for CI-barcodes to a unique number and
 *        save positions for CI-barocdes, AB, treatment barcode within the
 *        file of barcodes (position among all varying barcodes with [NNN...] pattern)
 */
void generateBarcodeDicts(std::string barcodeFile, std::string barcodeIndices, 
                          NBarcodeInformation& barcodeIdData, 
                          std::vector<std::string>& proteinDict, const int& protIdx, 
                          std::vector<std::string>* treatmentDict = nullptr, const int& treatmentIdx = 0);

/**
 * @brief A class to handle the processing of the demultiplexed data. 
 * This involves:
 *  - Mapping the barcodes to unique cell IDs, ABs, treatments
 *  - correcting UMIs, that are within a certain edit-distance
 *  - count finally all UMIs per single-cell AB combination
 *  - correct for errors in the data (same AB-scID has different UMIs, different treatments, etc. )
 */
class BarcodeProcessingHandler
{


    public:

        BarcodeProcessingHandler(NBarcodeInformation barcodeIdData) : varyingBarcodesPos(barcodeIdData){}

        void parseFile(const std::string fileName, const int& thread);

        //correct the umis, so that UnprocessedDemultiplexedData holds the same UMI for UMIs within a certain mismatch range
        //concurrently fills AbData, while iterating over AbSc lines, count the Ab occurence if this UMI has not been seen before
        void processBarcodeMapping(const int& umiMismatches, const int& thread);

        void writeStats(std::string output);
        void writeUmiCorrectedData(const std::string& output);
        inline void addTreatmentData(std::unordered_map<std::string, std::string > map)
        {
            rawData.setTreatmentDict(map);
        }
        inline void addProteinData(std::unordered_map<std::string, std::string > map)
        {
            rawData.setProteinDict(map);
        }

        void extended_umi_quality_check(const int& thread, const std::string& output);

    private:

        //parse the file, store each line in UnprocessedDemultiplexedData structure (ABs, treatment is already stored as a name,
        // single cells are defined by a dot seperated list of indices)
        void addFastqReadToUmiData(const std::string& line, const int& elements);
        void parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads);
        void correctUmis(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<scAbCount>& abDataTmp, 
                         const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected, 
                         const std::vector<dataLinePtr>& dataLinesToDelete);
        bool checkDataLineValidityDueToUmiBackground(const dataLinePtr& line, const std::vector<dataLinePtr>& dataLinesToDelete);

        void removeFalseSingleCellsFromUmis(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, int& currentUmisChecked, std::vector<dataLinePtr>& dataLinesToDelete);
        void correctUmisWithStats(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<scAbCount>& abDataTmp, 
                         const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected);

        void umiQualityCheck(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, umiQuality& qualTmp, int& currentUmisChecked);
        void umiQualityCheckExtended(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, int& currentUmisChecked, 
                                     std::vector<std::pair<unsigned long long, unsigned long long>>& uniqueUmiToDiffAbSc, 
                                     std::vector<umiQualityExtended>& extendedQuality, const std::string& output);

        //get positions of all barcodes in the lines of demultiplexed data
        void getBarcodePositions(const std::string& line, int& barcodeElements);

        //map all the barcodes of CI to a unique 'number' string as SingleCellIdx
        std::string generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes);

        //data structure storing lines with: UMI, AB_id, SingleCell_id
        //this is the raw data not UMI corrected
        UnprocessedDemultiplexedData rawData;

        //a final data structure without summed Ab-counts storing the corrected dataLines, after removal of erroneous 
        //dataLines and UMI corrections
        std::vector<dataLinePtr> umiData;

        //holding the counts per unique UMI (just for quality checks)
        std::vector<umiCount> newUmiData;
        //final data structures storing scID, AB-name, treatment-name and AB-count
        std::vector<scAbCount> abData;

        //statistics of the whole process
        StatsUmi stats;    
        umiQuality qual;  

        std::mutex lock;  
        
        //stores all the indices of variable barcodes in the barcode file (among all barcodes of [NNN...] pattern)
        NBarcodeInformation varyingBarcodesPos;
        // variables to read the data from each fastq line
        std::vector<int> fastqReadBarcodeIdx; // ids of the CI barcode within the whole pattern of all barcodes
        int abIdx = INT_MAX;
        int umiIdx = INT_MAX;
        int treatmentIdx = INT_MAX;
        int umiLength = 0;
};
