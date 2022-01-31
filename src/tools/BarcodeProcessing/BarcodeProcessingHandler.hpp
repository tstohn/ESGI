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
#include "helper.hpp"

/**
 * @brief Structure storing a vector with a mapping of the barcode-sequence to a unique ID
 *        for this sequence (this is done seperately for each barcoding round)
 * 
 *        It also stores also the indices for Ci-barcodes, AB, treatment.
 *        These indices indicate the position within the file of all barcodes. So only the position
 *        among all barcodes with the [NNN...] pattern - this data is used to then get the position in the
 *        full dataLine of Demultiplexed file.
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

//data type to represent final processed AB counts per single cell
struct scAbCount
{
    const char* abName;
    const char* treatment;
    
    const char* scID;
    int abCount = 0;
}; 

//data type representing counts per unique UMI in final processed data (without collapsed UMIs)
struct umiCount
{
    const char* umi;
    const char* abName;
    const char* treatment;
    
    const char* scID;
    int abCount = 0;
}; 

//some information about the read/ UMI quality (how many reads removed, how many Mismatches, etc.)
struct ProcessingLog
{
    //number of mismatches within reads for same AB/SC
    unsigned long long umiMM = 0; //UMIs with mismatches that were converted into another one 
                                    //(after reads were removed for 1.) non unique UMIs 2.) no class mapped to them )
    
    //number of removed reads for ...
    unsigned long long removedUmiReads = 0; // removed reads bcs. their UMI was not unique and was not present in >= 90% of the reads
    unsigned long long removedClassReads = 0; // removed reads bcs. the single cell had no class mapped (after non unique UMIs were already removed)
    
    //additional information about classes 
    unsigned long long removedSCDueToNoClass = 0; //equivalent for the reads with no class, but counting only unique cells that r removed
    uint removedClasses = 0; //number of mappings [SingleCell => Class] that were removed bcs.for a single cell all reads with guide did
                             // not represent to >=90% one unique Class

    //total number of raw reads, without the removed reads due to processing (the removed reads are removed from this full set of reads)
    unsigned long long totalReads = 0;
    unsigned long long totalGuideReads = 0;
    unsigned long long totalAbReads = 0;
};

/**
 * @brief Class storing all the final values (after removing erroneous reads with non unique UMIs, correcting MIsmatches in UMIs)
 *        the abData = final AB count data (represents a singleCell * AB matrix)
 *        the umiData = final UMI count data (same data as above but not UMI collapsed)
 *        the logData = basic values for number of removed reads due to mismatches, non-unique UMIs, etc.
 * 
 *        All the values are thread safe.
 */
class Results
{
    public:
        //getter functions (assume we only get those value AFTER they are completely filled -> no locking at this point)
        const std::vector<scAbCount> get_ab_data() const
        {
            return(abData);
        }
        const std::vector<umiCount> get_umi_data() const
        {
            return(umiData);
        }
        const ProcessingLog get_log_data() const
        {
            return(logData);
        }

        //setter functions, all locked for manipulation by threads
        void add_ab_count(const scAbCount& abCount)
        {
            abLock.lock();
            abData.push_back(abCount);
            abLock.unlock();
        }
        void add_umi_count(const umiCount& umiCount)
        {
            umiLock.lock();
            umiData.push_back(umiCount);
            umiLock.unlock();
        }
        void add_umi_mismatches(const unsigned long long& mm)
        {
            logLock.lock();
            ullong_save_add(logData.umiMM, mm);
            logLock.unlock();
        }
        void add_removed_reads_umi(const unsigned long long& reads)
        {
            logLock.lock();
            ullong_save_add(logData.removedUmiReads, reads);
            logLock.unlock();
        }
        void add_removed_reads_class(const unsigned long long& reads)
        {
            logLock.lock();
            ullong_save_add(logData.removedClassReads, reads);
            logLock.unlock();
        }
        void add_removed_class_for_single_cell()
        {
            logLock.lock();
            ++logData.removedClasses;
            logLock.unlock();
        }
        void set_total_reads(const unsigned long long& totalReads)
        {
            logLock.lock();
            logData.totalReads = totalReads;
            logLock.unlock();
        }
        void set_total_guide_reads(const unsigned long long& totalGuideReadsTmp)
        {
            logLock.lock();
            logData.totalGuideReads = totalGuideReadsTmp;
            logLock.unlock();
        }
        void set_total_ab_reads(const unsigned long long& totalAbReadsTmp)
        {
            logLock.lock();
            logData.totalAbReads = totalAbReadsTmp;
            logLock.unlock();
        }

    private:
        //holding the counts per unique UMI (just for quality checks)
        std::vector<umiCount> umiData;
        //final data structures storing scID, AB-name, treatment-name and AB-count
        std::vector<scAbCount> abData;
        //statistics of the whole process
        ProcessingLog logData;

        std::mutex umiLock;
        std::mutex abLock;
        std::mutex logLock;
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
 * 
 *  Overall workflow is:
 *  1.) parse lines (we cna parse normal AB counts and guide data for same single cells)
 *  2.) write temporary lines into for each umi (dict mapping umi to all their reads)
 *  3.) after collapsing UMI we write final liens into a dict mapping counts for AB+SC
 */
class BarcodeProcessingHandler
{


    public:

        BarcodeProcessingHandler(NBarcodeInformation barcodeIdData) : varyingBarcodesPos(barcodeIdData){}

        void parseFile(const std::string fileName, const int& thread);

        //counts AB and UMIs per single cell, data is stored in result (also saves basic information about processing
        //like removed reads, mismatched UMIs, etc.)
        //basic steps: 
        //1.) retain only reads for UMI has for more than 90% same Ab-Sc-Treatment
        //2.) retain only reads for SC where reads have more than 90% same treatment
        //3.) collapse same UMI reads with correcting mismatches in UMI for default 2 mismatches
        void processBarcodeMapping(const int& umiMismatches, const int& thread);

        void writeLog(std::string output);
        void writeAbCountsPerSc(const std::string& output);

        inline void addTreatmentData(std::unordered_map<std::string, std::string > map)
        {
            rawData.setTreatmentDict(map);
        }
        inline void addProteinData(std::unordered_map<std::string, std::string > map)
        {
            rawData.setProteinDict(map);
        }
        inline void addClassData(std::unordered_map<std::string, std::string > map)
        {
            rawData.setClassDict(map);
        }
        const UnprocessedDemultiplexedData getRawData() const
        {
            return rawData;
        }
        const std::vector<int> getCIBarcodeIdx() const
        {
            return(fastqReadBarcodeIdx);
        }
        const int getAbIdx() const
        {
            return(abIdx);
        }
        const int getUmiIdx() const
        {
            return(umiIdx);
        }
                const int getTreatmentIdx() const
        {
            return(treatmentIdx);
        }


    private:

        //parse the file, store each line in UnprocessedDemultiplexedData structure (ABs, treatment is already stored as a name,
        // single cells are defined by a dot seperated list of indices)
        void add_line_to_temporary_data(const std::string& line, const int& elements,
                                        std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>>& scClasseCountDict,
                                        unsigned long long& abReadCount, unsigned long long& guideReadCount);
        void parseBarcodeLines(std::istream* instream, const unsigned long long& totalReads, unsigned long long& currentReads,
                               std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>>& scClasseCountDict);
        
        //check if a read is in 'dataLinesToDelete' (not-unique UMI for this read)
        bool checkIfLineIsDeleted(const dataLinePtr& line, const std::vector<dataLinePtr>& dataLinesToDelete);
        //write all reads which come from a not-unique UMI into 'dataLinesToDelete' vector
        //these lines are later not considered when the AB-count per single cell is calculated
        void markReadsWithNoUniqueUmi(const std::vector<umiDataLinePtr>& uniqueUmis,
                                      std::atomic<unsigned long long>& count,
                                      const unsigned long long& totalCount);
        void markReadsWithNoUniqueTreatment(const std::vector<dataLinePtr>& uniqueSc,
                                                              std::vector<dataLinePtr>& dataLinesToDelete, 
                                                              std::atomic<unsigned long long>& count,
                                                              const unsigned long long& totalCount);
        //used within 'count_abs_per_single_cell' to get counts per UMI for reads of one AB SC combination
        void count_umi_occurence(std::vector<int>& positionsOfSameUmi, 
                                                   umiCount& umiLineTmp,
                                                   const std::vector<dataLinePtr>& allScAbCounts,
                                                   const int& umiMismatches,
                                                   const int& lastIdx);
        //count the ABs per single cell (iterating over reads for a AB-SC combination and summing them)
        //(not considering not-unique UMIs in 'dataLinesToDelete' vector, collapsing UMIs,etc.)
        void count_abs_per_single_cell(const int& umiMismatches, const std::vector<dataLinePtr>& uniqueAbSc,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount,
                                                        std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
                    CharHash, CharPtrComparator>>& umiMap);

        //get positions of all barcodes in the lines of demultiplexed data
        void getBarcodePositions(const std::string& line, int& barcodeElements);

        //map all the barcodes of CI to a unique 'number' string as SingleCellIdx
        std::string generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes);

        //functions processing the class labels for single cells (obtained by guide reads)
        void generate_unique_sc_to_class_dict(const std::unordered_map< const char*, 
                                              std::unordered_map< const char*, unsigned long long>>& scClasseCountDict);

        void combine_ab_and_guide_data(std::vector<dataLine>& abDataLines,
                                       const std::unordered_map< const char*, const char*>& scClassMap);

        //data structure storing lines with: UMI, AB_id, SingleCell_id
        //this is the raw data not UMI corrected
        UnprocessedDemultiplexedData rawData;
        // the final data: ABCounts, UMICounts, and a processingLog containing basic values (removed reads, etc.)
        Results result;

        std::mutex statusUpdateLock;  
        std::mutex writeToRawDataLock; //while processing reads of same UMI, we write UMI collapsed reads into the dict
        //for reads of same AB/SC and have to lock writing

        // DATA STRUCTURES FOR PARSING DEMULTIPLEXED DATA
        //stores all the indices of variable barcodes in the barcode file (among all barcodes of [NNN...] pattern)
        NBarcodeInformation varyingBarcodesPos;
        // variables to read the data from each fastq line
        std::vector<int> fastqReadBarcodeIdx; // ids of the CI barcode within the whole pattern of all barcodes
        int abIdx = INT_MAX;
        int umiIdx = INT_MAX;
        int treatmentIdx = INT_MAX;
        int umiLength = 0;
};
