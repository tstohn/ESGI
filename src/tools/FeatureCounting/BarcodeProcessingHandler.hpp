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
#include <regex>

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
 *        It stores the indices (0-indexed) that contain information of the whole tsv file
 *        (an index is the position in the barcode tsv file)
 */
struct BarcodeInformation
{
    // mapping ColumnIdx -> map(barcodes -> number)
    //these dicts are in order of scBarcodeIndices!!!!
    std::vector<std::unordered_map<std::string, int>> barcodeIdMaps; //used to generate a SC-index
    
    //different indices, they r the index of the NNN-barcodes only (and therefore differe from the index of all barcodes incl. constant ones)
    std::vector<int> scBarcodeIndices; //CI barcode indices
    int treatmentIdx = -1; // treatment index
    unsigned int featureIdx = 1; //AB index
    
    std::vector<int> umiIdx; //UMI index
    unsigned int umiMismatches;
    unsigned int umiLength = 0; //the sum of lengths for all UMIs
};

//data type to represent final processed AB counts per single cell
struct scAbCount
{
    const char* abName;
    const char* treatment;
    const char* className;

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

struct umiDist
{
    std::unordered_map<std::string, int> abs;
    std::vector<std::unordered_map<int, unsigned long>> dist;
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
    unsigned int removedClasses = 0; //number of mappings [SingleCell => Class] that were removed bcs.for a single cell all reads with guide did
                             // not represent to >=90% one unique Class

    //total number of raw reads, without the removed reads due to processing (the removed reads are removed from this full set of reads)
    unsigned long long totalReads = 0;
    unsigned long long totalGuideReads = 0;
    unsigned long long totalAbReads = 0;
};

//functions to open gzipped/ non-zipped files containing the demultiplexed reads
inline std::istream* openFile(const std::string& filename,
                       std::ifstream& fileStream,
                       boost::iostreams::filtering_streambuf<boost::iostreams::input>& inbuf,
                       bool isGz)
{
    fileStream.close();  // In case it's already open
    fileStream.clear();  // Clear any EOF flags
    fileStream.open(filename, std::ios_base::in | std::ios_base::binary);
    if (!fileStream.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return nullptr;
    }

    if (isGz) 
    {
        inbuf.reset();  // Clear previous filters
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(fileStream);
        return new std::istream(&inbuf);
    } else {
        return &fileStream;
    }
}
inline bool isGzipped(const std::string& filename) 
{
    return boost::algorithm::iends_with(filename, ".gz");
}

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
        const umiDist get_umi_stats() const
        {
            return(umiStat);
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
        void add_umi_stats(const umiCount& umiCount)
        {
            statLock.lock();
            std::unordered_map<std::string, int>::iterator iter = umiStat.abs.find(umiCount.abName);
            if(iter == umiStat.abs.end())
            {
                iter = umiStat.abs.insert(std::make_pair<std::string,int>(umiCount.abName, umiStat.abs.size())).first;
                std::unordered_map<int, unsigned long> umap;
                umiStat.dist.push_back(umap);
            }

            std::unordered_map<int, unsigned long>* AbDict = &umiStat.dist.at(iter->second);
            std::unordered_map<int, unsigned long>::iterator iterStat = AbDict->find(umiCount.abCount);

            if(iterStat == AbDict->end())
            {
                (*AbDict)[umiCount.abCount] = 1;
            }
            else
            {
                (*AbDict)[iterStat->first] += 1;
            }

            statLock.unlock();
        }

    private:
        //holding the counts per unique UMI (just for quality checks)
        std::vector<umiCount> umiData;
        //final data structures storing scID, AB-name, treatment-name and AB-count
        std::vector<scAbCount> abData;
        //statistics of the whole process
        ProcessingLog logData;

        //statistics of UMIs. maps the amplifcation of a umi to the occurence in data
        umiDist umiStat;

        std::mutex umiLock;
        std::mutex abLock;
        std::mutex logLock;
        std::mutex statLock;
};

/**
 * @brief Map barcode sequences for CI-barcodes to a unique number and
 *        save positions for CI-barocdes, AB, treatment barcode within the
 *        file of barcodes (position among all varying barcodes with [NNN...] pattern)
 */
void generateBarcodeDicts(const std::string& headerLine, const std::string& barcodeDir, std::string barcodeIndices, 
                          BarcodeInformation& barcodeIdData, 
                          std::vector<std::string>& proteinNamelist, bool parseAbBarcodes, const int& featureIdx, 
                          std::vector<std::string>* treatmentDict = nullptr, const int& treatmentIdx = -1,
                          std::string umiIdx = "", int umiMismatches = 1);

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

        BarcodeProcessingHandler(BarcodeInformation barcodeInformationInput) : barcodeInformation(barcodeInformationInput){}

        void parse_barcode_file(std::string& inFile);

        //counts AB and UMIs per single cell, data is stored in result (also saves basic information about processing
        //like removed reads, mismatched UMIs, etc.)
        //basic steps: 
        //1.) retain only reads for UMI has for more than 90% same Ab-Sc-Treatment
        //2.) retain only reads for SC where reads have more than 90% same treatment
        //3.) collapse same UMI reads with correcting mismatches in UMI for default 2 mismatches
        void processBarcodeMapping(const int& thread);

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
            return(barcodeInformation.scBarcodeIndices);
        }
        int getAbIdx() const
        {
            return(barcodeInformation.featureIdx);
        }
        const std::vector<int> getUmiIdx() const
        {
            return(barcodeInformation.umiIdx);
        }
        int getTreatmentIdx() const
        {
            return(barcodeInformation.treatmentIdx);
        }
        void setUmiFilterThreshold(double threshold)
        {
            umiFilterThreshold = threshold;
        }
        void setScClassConstaint(bool scMustHaveClassTmp)
        {
            scMustHaveClass = scMustHaveClassTmp;
        }
        void setumiRemoval(bool umiRemovalTmp)
        {
            umiRemoval = umiRemovalTmp;
        }
        void setSingleCellIdStyle(bool scIdStringTmp)
        {
            scIdString = scIdStringTmp;
        }

    private:

        //parse the file, store each line in UnprocessedDemultiplexedData structure (ABs, treatment is already stored as a name,
        // single cells are defined by a dot seperated list of indices)
        void add_line_to_temporary_data(const std::string& line, const size_t& elements,
                                        unsigned long long& readCount);
        void parseBarcodeLines(std::string& inFile, const unsigned long long& totalReads, unsigned long long& currentReads);
        
        //check if a read is in 'dataLinesToDelete' (not-unique UMI for this read)
        bool checkIfLineIsDeleted(const dataLinePtr& line, const std::vector<dataLinePtr>& dataLinesToDelete);
        //stores a real unique read in a dict for the corresponding AB-SC (only read with UMI presence > 90 considered)
        //reads r collapsed
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
                                                   const int& lastIdx);
        //count the ABs per single cell (iterating over reads for a AB-SC combination and summing them, this is already a sparse vector)
        //reads of same UMI are collapsed before
        void count_abs_per_single_cell(const std::vector<dataLinePtr>& uniqueAbSc,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount,
                                                        std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
                                                        CharHash, CharPtrComparator>>& umiMap);

        //get positions of all barcodes in the lines of demultiplexed data
        void getBarcodePositions(const std::string& line, int& barcodeElements);

        //map all the barcodes of CI to a unique 'number' string as SingleCellIdx
        std::string generateSingleCellIndexFromBarcodes(const std::vector<std::string>& ciBarcodes);

        void combine_ab_and_guide_data(std::vector<dataLine>& abDataLines,
                                       const std::unordered_map< const char*, const char*>& scClassMap);

        //data structure storing lines with: UMI, AB_id, SingleCell_id
        //this is the raw data not UMI corrected
        UnprocessedDemultiplexedData rawData;
        // the final data: ABCounts, UMICounts, and a processingLog containing basic values (removed reads, etc.)
        Results result;
        std::unordered_map< const char*, unsigned long long> guideCountPerSC;

        std::mutex statusUpdateLock;  
        std::mutex writeToRawDataLock; //while processing reads of same UMI, we write UMI collapsed reads into the dict
        //for reads of same AB/SC and have to lock writing

        // DATA STRUCTURES FOR PARSING DEMULTIPLEXED DATA
        //stores all the indices of variable barcodes in the barcode file (among all barcodes of [NNN...] pattern)
        BarcodeInformation barcodeInformation;

        double umiFilterThreshold = 0.0;
        bool scMustHaveClass = true;
        bool umiRemoval = true;
        bool scIdString = false;
};
